#!/usr/bin/env python3
"""
vel_loop_id_v3.py -- identification du retard de la boucle de vitesse ODrive.

CE QUE MESURE CE SCRIPT, ET POURQUOI
====================================
Toute la synthese LQR/TVLQR repose sur

    theta_ddot = u

c'est-a-dire sur l'hypothese que la boucle de vitesse du variateur est
INFINIMENT RAPIDE : la consigne envoyee est la vitesse obtenue. C'est cette
hypothese qui autorise control_loop a alimenter le retour d'etat avec
l'integrale de u (USE_MEASURED_THETA = False) plutot qu'avec l'encodeur.

Les journaux disent ou cette hypothese tient et ou elle cesse de tenir :

    campagne          RMS |vel_cmd - vel_enc|      max
    T1 / T2               0.15 - 0.40 rad/s      0.5 - 1.0
    T3 (looping)          1.50 - 2.72 rad/s      4.6 - 7.7

En T3 l'ecart atteint 20 % de la consigne. Or theta_int est l'INTEGRALE de cet
ecart : sur les 2.5 a 3.1 s d'un plan, un biais moyen de 1 rad/s laisse
plusieurs radians d'erreur sur theta, a comparer a un |theta_r| planifie du
meme ordre. Le terme K1*(theta - theta_ref) du TVLQR regule donc sur une
grandeur qui a cesse de decrire le cerceau, et il le fait a pleine autorite.

Deux facons d'en sortir. Lire l'encodeur a chaque tour coute du temps USB et
referme une boucle parasite autour de la dynamique du variateur, non modelisee
-- c'est exactement ce que l'en-tete de control_loop redoutait, a juste titre.
L'autre est de garder l'integration en boucle ouverte mais de lui faire porter
le retard, en filtrant la consigne par un premier ordre :

    theta_dot_hat[k] = theta_dot_hat[k-1]
                     + dt/(TAU_VEL_S + dt) * (theta_dot_cmd[k] - theta_dot_hat[k-1])

Cout : zero aller-retour USB, zero bruit d'encodeur. Encore faut-il connaitre
TAU_VEL_S, et c'est l'objet de ce fichier.

CE QUI EST IDENTIFIE, EXACTEMENT
================================
La reponse indicielle mesuree ICI, depuis le Raspberry Pi, contient
inseparablement :
  - la dynamique du regulateur de vitesse de l'ODrive,
  - la constante electrique et mecanique du groupe moteur + cerceau,
  - le retard d'ecriture USB de la consigne,
  - le retard de lecture USB de vel_estimate.
C'est une bonne chose : c'est precisement le retard que la loi de commande
subit, donc celui qu'elle doit modeliser. Ce n'est PAS la bande passante
interne du variateur, et ce chiffre ne doit pas etre presente comme telle.

DEUX MODELES SONT AJUSTES
=========================
  1. premier ordre pur      v(t) = V (1 - exp(-t/tau_eq))
     C'est la structure du filtre implante dans control_loop, donc tau_eq est
     LA valeur a reporter dans TAU_VEL_S. Rien d'autre n'est directement
     utilisable par le code.

  2. premier ordre + retard pur (FOPDT)
     v(t) = V (1 - exp(-(t-L)/tau))  pour t >= L, 0 sinon
     Il separe ce qui est du transport (L, essentiellement l'USB et la
     quantification de la consigne a 50 Hz) de ce qui est de la dynamique
     (tau). Cette decomposition ne sert pas au code mais elle est ce qu'il
     faut ecrire dans un memoire : elle dit ou va le temps.
     [Ljung, "System Identification: Theory for the User", 2e ed., 1999,
      ch. 6 ; Astrom & Hagglund, "PID Controllers", 2e ed., 1995, ch. 2.]

L'ajustement porte sur PLUSIEURS echelons, dans LES DEUX SENS, et la dispersion
entre echelons est reportee : un modele du premier ordre qui ne rendrait pas
compte du frottement sec donnerait des tau systematiquement differents selon le
sens, et il faut le voir plutot que le moyenner en silence.

PRECAUTIONS
===========
Le cerceau TOURNE pendant cet essai. Retirer la balle. L'amplitude par defaut
est volontairement modeste (5 rad/s, soit 0.8 tr/s) : la constante de temps
d'un systeme lineaire n'en depend pas, et une amplitude plus faible reduit le
risque materiel. Si tau varie nettement avec l'amplitude, le systeme n'est pas
lineaire dans cette plage et c'est un resultat en soi.

Usage
-----
    T1_HOOP=outer python3 vel_loop_id_v3.py
    T1_HOOP=outer python3 vel_loop_id_v3.py --amp 10 --reps 5
    T1_HOOP=outer python3 vel_loop_id_v3.py --amp 5 --amp2 15   # test de linearite
"""

import argparse
import json
import math
import os
import time

import numpy as np

from common import bench_common as bc


# =========================================================================== #
#  Acquisition d'un echelon
# =========================================================================== #

def record_step(ax, link, v_from, v_to, hold_s):
    """Envoie un echelon de consigne et journalise vel_estimate au plus vite.

    Aucune cadence n'est imposee : la boucle lit aussi vite que l'USB le
    permet (typiquement 300-600 Hz), ce qui est necessaire pour resoudre une
    constante de temps attendue autour de 15 ms. Les instants sont ceux du
    RETOUR de la lecture, donc le retard de lecture est inclus -- voir
    l'en-tete : c'est voulu.
    """
    t, v = [], []
    t0 = time.perf_counter()
    link.send_velocity(v_to)
    while True:
        now = time.perf_counter() - t0
        if now >= hold_s:
            break
        vel = ax.vel_estimate * 2.0 * math.pi
        t.append(time.perf_counter() - t0)
        v.append(vel)
    return np.asarray(t), np.asarray(v)


# =========================================================================== #
#  Ajustements
# =========================================================================== #

def _normalise(t, v, v_from):
    """Ramene l'echelon a une transition 0 -> 1.

    La valeur finale est prise comme la MEDIANE du dernier tiers : la moyenne
    serait tiree par le transitoire si le maintien est court, et la mediane
    est insensible aux lectures aberrantes.
    """
    n = t.size
    v_end = float(np.median(v[int(0.67 * n):]))
    span = v_end - v_from
    if abs(span) < 1e-6:
        return None, None, v_end
    return t, (v - v_from) / span, v_end


def fit_first_order(t, y):
    """v(t) = 1 - exp(-t/tau). Retourne tau, ou None si l'ajustement echoue.

    Ajustement non lineaire par moindres carres. L'initialisation vient du
    temps de montee a 63.2 %, ce qui est deja l'estimateur graphique
    classique : le raffinement ne fait que le rendre insensible au bruit d'un
    echantillon unique.
    """
    from scipy.optimize import curve_fit

    def model(tt, tau, gain):
        return gain * (1.0 - np.exp(-np.clip(tt / max(tau, 1e-6), 0, 50)))

    i = np.argmax(y >= 0.632)
    tau0 = max(float(t[i]), 1e-3) if y.max() >= 0.632 else 0.02
    try:
        popt, _ = curve_fit(model, t, y, p0=[tau0, 1.0],
                            bounds=([1e-4, 0.5], [1.0, 1.5]), maxfev=20000)
    except Exception:
        return None, None
    resid = y - model(t, *popt)
    return float(popt[0]), float(np.sqrt(np.mean(resid ** 2)))


def fit_fopdt(t, y):
    """v(t) = 1 - exp(-(t-L)/tau) pour t >= L. Retourne (L, tau, rms)."""
    from scipy.optimize import curve_fit

    def model(tt, L, tau, gain):
        z = np.clip((tt - L) / max(tau, 1e-6), 0.0, 50.0)
        return gain * (1.0 - np.exp(-z))

    try:
        popt, _ = curve_fit(model, t, y, p0=[0.005, 0.015, 1.0],
                            bounds=([0.0, 1e-4, 0.5], [0.2, 1.0, 1.5]),
                            maxfev=20000)
    except Exception:
        return None, None, None
    resid = y - model(t, *popt)
    return (float(popt[0]), float(popt[1]),
            float(np.sqrt(np.mean(resid ** 2))))


# =========================================================================== #
#  Campagne
# =========================================================================== #

def run_campaign(ax, link, amp, reps, hold_s):
    """Sequence 0 -> +A -> 0 -> -A -> 0, repetee. Retourne la liste des
    transitions exploitables."""
    targets = []
    for _ in range(reps):
        targets += [amp, 0.0, -amp, 0.0]

    results = []
    v_from = 0.0
    for k, v_to in enumerate(targets):
        t, v = record_step(ax, link, v_from, v_to, hold_s)
        tt, y, v_end = _normalise(t, v, v_from)
        if tt is not None:
            tau, rms1 = fit_first_order(tt, y)
            L, tau2, rms2 = fit_fopdt(tt, y)
            results.append({
                "from": v_from, "to": v_to, "v_end": v_end,
                "n": int(t.size), "rate_hz": float(t.size / hold_s),
                "tau_eq": tau, "rms_1st": rms1,
                "L": L, "tau_fopdt": tau2, "rms_fopdt": rms2,
            })
            print("    {:2d}  {:+6.2f} -> {:+6.2f} rad/s  "
                  "({:5.0f} Hz, {:4d} pts)   tau_eq = {}"
                  .format(k + 1, v_from, v_to, t.size / hold_s, t.size,
                          "  ---  " if tau is None
                          else "{:5.1f} ms".format(1e3 * tau)))
        v_from = v_to
    return results


def summarise(results, label):
    """Statistiques robustes et separation par sens de l'echelon."""
    ok = [r for r in results if r["tau_eq"] is not None]
    if not ok:
        print("\n[!!] aucun echelon exploitable pour {}".format(label))
        return None

    tau = np.array([r["tau_eq"] for r in ok])
    up = np.array([r["tau_eq"] for r in ok if r["to"] > r["from"]])
    dn = np.array([r["tau_eq"] for r in ok if r["to"] < r["from"]])
    Lv = np.array([r["L"] for r in ok if r["L"] is not None])
    tf = np.array([r["tau_fopdt"] for r in ok if r["tau_fopdt"] is not None])

    print("\n  --- {} ({} echelons exploitables) ---".format(label, len(ok)))
    print("    premier ordre pur")
    print("      tau_eq        : {:5.1f} ms   (mediane)".format(1e3 * np.median(tau)))
    print("      dispersion    : {:5.1f} ms   (ecart-type)".format(1e3 * tau.std()))
    print("      etendue       : {:5.1f} - {:5.1f} ms"
          .format(1e3 * tau.min(), 1e3 * tau.max()))
    if up.size and dn.size:
        print("      montee / desc.: {:5.1f} / {:5.1f} ms"
              .format(1e3 * np.median(up), 1e3 * np.median(dn)))
        ecart = abs(np.median(up) - np.median(dn))
        if ecart > 0.3 * np.median(tau):
            print("      [!] {:.0f} % d'ecart entre les deux sens : le "
                  "frottement sec".format(100 * ecart / np.median(tau)))
            print("          n'est pas negligeable, le premier ordre est une "
                  "approximation.")
    print("      residu moyen  : {:.4f} (sur une transition normalisee a 1)"
          .format(float(np.mean([r["rms_1st"] for r in ok]))))

    if Lv.size and tf.size:
        print("    premier ordre + retard pur")
        print("      L (transport) : {:5.1f} ms".format(1e3 * np.median(Lv)))
        print("      tau (dynamique): {:5.1f} ms".format(1e3 * np.median(tf)))
        print("      L + tau       : {:5.1f} ms   (a comparer a tau_eq)"
              .format(1e3 * (np.median(Lv) + np.median(tf))))

    return float(np.median(tau))


# =========================================================================== #
#  Programme principal
# =========================================================================== #

def main():
    from tasks import t3_config as cfg

    p = argparse.ArgumentParser(
        description="Identification du retard de la boucle de vitesse ODrive")
    p.add_argument("--amp", type=float, default=5.0,
                   help="amplitude de l'echelon [rad/s]")
    p.add_argument("--amp2", type=float, default=None,
                   help="seconde amplitude, pour tester la linearite")
    p.add_argument("--reps", type=int, default=4,
                   help="nombre de cycles (4 transitions par cycle)")
    p.add_argument("--hold", type=float, default=0.40,
                   help="duree de maintien apres chaque echelon [s]")
    p.add_argument("--outdir", default="data/id_vel")
    args = p.parse_args()

    if abs(args.amp) > cfg.THETA_DOT_MAX:
        raise SystemExit("amplitude au-dela de THETA_DOT_MAX = {:.0f} rad/s"
                         .format(cfg.THETA_DOT_MAX))

    print("=" * 70)
    print("  Identification de la boucle de vitesse -- cerceau {}"
          .format(cfg.HOOP))
    print("=" * 70)
    print("  amplitude    : {:+.1f} rad/s{}"
          .format(args.amp, "" if args.amp2 is None
                  else " puis {:+.1f}".format(args.amp2)))
    print("  cycles       : {} ({} transitions)".format(args.reps, 4 * args.reps))
    print("  maintien     : {:.2f} s par echelon".format(args.hold))
    print("\n  [!] LE CERCEAU VA TOURNER. Retirer la balle et degager le banc.")
    input("      ENTREE quand c'est fait, Ctrl-C pour abandonner ...")

    odrv0, ax = bc.connect_and_prepare(cfg)
    link = bc.ODriveLink(ax, cfg)

    tau_a = tau_b = None
    try:
        print("\n[id] campagne a {:+.1f} rad/s".format(args.amp))
        res_a = run_campaign(ax, link, args.amp, args.reps, args.hold)
        tau_a = summarise(res_a, "amplitude {:+.1f} rad/s".format(args.amp))

        res_b = []
        if args.amp2 is not None:
            print("\n[id] campagne a {:+.1f} rad/s".format(args.amp2))
            res_b = run_campaign(ax, link, args.amp2, args.reps, args.hold)
            tau_b = summarise(res_b,
                              "amplitude {:+.1f} rad/s".format(args.amp2))
    finally:
        bc.ramp_down(ax, 0.0, cfg)
        print("\n[odrv] cerceau arrete")

    if tau_a is None:
        raise SystemExit("identification echouee")

    if tau_b is not None:
        ecart = abs(tau_b - tau_a) / max(tau_a, 1e-9)
        print("\n  --- linearite ---")
        print("    tau_eq({:+.1f}) = {:.1f} ms, tau_eq({:+.1f}) = {:.1f} ms, "
              "ecart {:.0f} %".format(args.amp, 1e3 * tau_a,
                                      args.amp2, 1e3 * tau_b, 100 * ecart))
        if ecart > 0.25:
            print("    [!] la constante de temps depend de l'amplitude : le")
            print("        premier ordre lineaire n'est valable que par "
                  "morceaux.")
            print("        Retenir la valeur mesuree a l'amplitude la plus")
            print("        proche de celle du plan T3.")

    tau = tau_a if tau_b is None else max(tau_a, tau_b)
    print("\n  --- a recopier dans t3_config_v3.py ---")
    print("    TAU_VEL_S = {:.4f}     # {:.1f} ms".format(tau, 1e3 * tau))
    print("\n  Rappel : TAU_VEL_S = 0.0 reproduit exactement le comportement")
    print("  T1/T2 anterieur. Ne le renseigner que pour T3, et le laisser a 0")
    print("  ailleurs tant que les campagnes T1/T2 ne sont pas refaites.")

    os.makedirs(args.outdir, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    path = os.path.join(args.outdir, "vel_loop_id_{}_{}.json"
                        .format(cfg.HOOP, stamp))
    with open(path, "w") as f:
        json.dump({"hoop": cfg.HOOP, "timestamp": stamp,
                   "amp": args.amp, "amp2": args.amp2,
                   "hold_s": args.hold, "reps": args.reps,
                   "tau_vel_s": tau,
                   "steps_amp": res_a, "steps_amp2": res_b}, f, indent=2)
    print("\n[out] {}".format(path))


if __name__ == "__main__":
    main()
