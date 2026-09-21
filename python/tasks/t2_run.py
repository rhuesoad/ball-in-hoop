#!/usr/bin/env python3
"""
t2_run_v3.py -- T2 : poursuite d'une reference variable par la balle.

Meme structure que t1_run_v3.py : toute la mecanique est dans
bench_common.py, ce fichier ne fournit que le generateur de reference. La
SEULE difference entre T1 et T2 est la nature de ce generateur.


POURQUOI UNE REFERENCE VARIABLE N'EST PAS UNE CONSIGNE VARIABLE
===============================================================

La loi u = -K(x - x_ref) n'a de sens que si le couple (x_ref, u_ref) est
une SOLUTION de la dynamique : il doit exister une trajectoire de
reference telle que

    xdot_r = A x_r + B u_r,     C x_r = psi_r(t)

La loi correcte est alors a deux degres de liberte

    u = u_r(t) - K (x - x_r(t))

ou u_r(t) PRODUIT le mouvement et K ne fait que ramener l'ecart a zero.
Un gain plus grand ne remplace jamais un u_r manquant.
[Franklin, Powell & Emami-Naeini, "Feedback Control of Dynamic Systems",
 ch. 7, gains N_x / N_u ; Astrom & Murray, "Feedback Systems", 2e ed.,
 ch. 7-8, structure feedforward + feedback.]

La version v1 de T2 posait x_ref = [0, 0, psi_ref(t), 0] et u_ff = 0.
Trois consequences, toutes mesurables :

  1. psi_dot_ref = 0 alors que la consigne bouge : le terme K4 agit en
     amortissement PUR CONTRE le mouvement demande. C'est la cause
     dominante de la perte d'amplitude.
  2. theta_ref = theta_dot_ref = 0 : le retour d'etat combat la rotation
     du cerceau qui est justement necessaire pour produire psi_r. A
     5 deg / 1 Hz le terme parasite K2*theta_dot est du meme ordre que le
     terme utile K3*(psi - psi_r).
  3. aucun feedforward : la boucle doit tout produire par l'erreur, donc
     l'erreur ne peut pas etre petite.

Ici x_r(t) a ses QUATRE composantes et u_ff(t) est calcule par inversion
du modele.

CE QUI SUBSISTE MALGRE LE FEEDFORWARD
-------------------------------------
Le feedforward annule l'erreur dans le cas NOMINAL. Toute erreur de modele
(A21, A22, B2, retard) laisse une erreur residuelle que le retour d'etat
seul ne peut pas annuler : un retour d'etat statique a un gain de boucle
FINI a toute frequence non nulle. L'annulation asymptotique exigerait que
la boucle contienne un modele du generateur de la reference, c'est-a-dire
une paire de poles a +/- j*omega (principe du modele interne, Francis &
Wonham, "The internal model principle of control theory", Automatica
12(5):457-465, 1976). Un integrateur classique ne suffit pas : son gain
n'est infini qu'en continu, pas a 1 Hz.

ADMISSIBILITE
-------------
Une consigne constante non nulle n'est pas admissible (theta_dot diverge,
voir t1_run_v3.py). Une sinusoide l'est, tant que l'effort demande au
cerceau reste sous THETA_DOT_MAX. Le script calcule cet effort EXACTEMENT
avant de lancer et refuse si la marge est insuffisante.

CHOIX DE L'ESTIMATEUR DE PSI_DOT
--------------------------------
    derivative  difference arriere + passe-bas du premier ordre.
    kalman      estimator_v3.DelayCompensatedKF, modele linearise psi = 0.
    ekf         ekf_v3.ExtendedKF, rappel en A21*sin(psi), jacobienne
                relinearisee a chaque pas. Identique au precedent tant que
                psi reste petit, ce qui est le cas en T2 : c'est justement
                ce qui permet de les comparer sans autre variable.

Usage
-----
    python3 t2_run_v3.py                                 # 5 deg, 1 Hz
    python3 t2_run_v3.py --amp 8 --freq 1.3              # pres de la resonance
    python3 t2_run_v3.py --no-feedforward                # comparaison v1
    python3 t2_run_v3.py --estimator ekf                 # KF etendu
    python3 t2_run_v3.py --ref step --amp 5              # non admissible
    python3 t2_run_v3.py --envelope                      # table, sans lancer
"""

import argparse
import math

import numpy as np

from common import bench_common as bc
from tasks import t2_config as cfg


def build_estimator(name, cfg):
    """Retourne (estimateur, n_delay) ou (None, 0) pour la derivee filtree."""
    if name == "derivative":
        return None, 0
    if name == "kalman":
        from common import estimator as est_mod
    elif name == "ekf":
        from common import ekf as est_mod
    else:
        raise SystemExit("estimateur inconnu : {}".format(name))
    return est_mod.build_estimator(cfg, cfg.DT_NOMINAL)


# =========================================================================== #
#  Generateurs de reference
# =========================================================================== #

def _fade(t):
    """Demi-cosinus, 0 -> 1 sur REF_FADE_S. Derivable, derivee nulle aux
    deux bouts."""
    if cfg.REF_FADE_S <= 0.0:
        return 1.0
    if t >= cfg.REF_FADE_S:
        return 1.0
    return 0.5 * (1.0 - math.cos(math.pi * t / cfg.REF_FADE_S))


def sine_reference(amp_deg, freq_hz, feedforward=True):
    """Sinusoide, trajectoire de reference COMPLETE par inversion du modele.

        psi_r      =  A sin(wt)
        psi_dot_r  =  A w cos(wt)
        psi_ddot_r = -A w^2 sin(wt)

        u_r = (psi_ddot_r - A21 psi_r - A22 psi_dot_r) / B2
            = P sin(wt) + Q cos(wt)

        theta_dot_r = (-P cos(wt) + Q sin(wt)) / w        (moyenne nulle)
        theta_r     = -(P sin(wt) + Q cos(wt)) / w^2      (moyenne nulle)

    Verification : d^2(theta_r)/dt^2 = P sin(wt) + Q cos(wt) = u_r. Les
    constantes d'integration sont choisies nulles pour que le cerceau
    oscille autour de sa position d'armement au lieu de deriver.

    PORTEE EXACTE DU FONDU
    ----------------------
    Le fondu g(t) multiplie x_r ET u_ff, mais ses derivees n'entrent pas
    dans l'inversion du modele. Le couple (g*x_r, g*u_r) n'est donc une
    solution EXACTE de la dynamique que pour g constant, c'est-a-dire
    apres REF_FADE_S. Pendant le fondu il manque les termes en gpoint et
    gseconde, d'ordre A*gpoint*w et A*gseconde par rapport a A*w^2 : pour
    REF_FADE_S = 2 s et f = 1 Hz, gpoint_max = pi/(2*T) = 0.79 s^-1 contre
    w = 6.28, soit une erreur relative de feedforward d'environ 25 % au
    milieu du fondu, qui decroit a zero a la fin.

    Ce n'est pas corrige, et volontairement : rendre le fondu exact
    imposerait d'integrer theta_r en ligne, donc d'accepter une derive
    lente de la position de reference sur toute la duree de l'essai -- un
    defaut permanent en echange d'un defaut transitoire. Le prix a payer
    est que les DEUX PREMIERES SECONDES d'un essai T2 ne sont pas
    exploitables pour chiffrer la poursuite. analyse_tracking() n'utilise
    que la seconde moitie de l'essai, donc la mesure publiee est propre ;
    c'est le trace qu'il faut lire en connaissant cette limite.
    """
    A = math.radians(amp_deg)
    w = 2.0 * math.pi * freq_hz
    P = -A * (w ** 2 + cfg.A21) / cfg.B2
    Q = -A * cfg.A22 * w / cfg.B2

    def ref(t):
        s, c = math.sin(w * t), math.cos(w * t)
        g = _fade(t)
        psi_r = A * s
        psi_dot_r = A * w * c
        if not feedforward:
            # Reproduit exactement la v1 : seule psi_ref est renseignee.
            return np.array([0.0, 0.0, g * psi_r, 0.0]), 0.0
        u_r = P * s + Q * c
        theta_dot_r = (-P * c + Q * s) / w
        theta_r = -(P * s + Q * c) / w ** 2
        return (np.array([g * theta_r, g * theta_dot_r,
                          g * psi_r, g * psi_dot_r]), g * u_r)
    return ref


def step_reference(amp_deg, feedforward=True):
    """Echelon. NON ADMISSIBLE : psi = cste != 0 impose theta_ddot constant,
    donc theta_dot diverge. On garde donc theta_r = theta_dot_r = 0 et
    u_ff = 0 (retour d'etat seul), et l'erreur statique observee EST le
    resultat de l'essai. Conserve uniquement comme point de comparaison
    avec T1."""
    A = math.radians(amp_deg)

    def ref(t):
        g = _fade(t)
        return np.array([0.0, 0.0, g * A, 0.0]), 0.0
    return ref


def triangle_reference(amp_deg, freq_hz, feedforward=True):
    """Triangle DEMARRANT A ZERO (la v1 partait a -A, soit un saut de
    consigne a la fermeture de boucle).

    psi_ddot_r est nulle par morceaux et impulsionnelle aux sommets : le
    feedforward exact n'existe pas. On utilise le feedforward partiel
    u_r = (-A21 psi_r - A22 psi_dot_r)/B2, et theta_r / theta_dot_r sont
    obtenus par integration numerique en ligne. La derive lente de cette
    integration est le prix a payer ; a n'utiliser que pour illustrer le
    comportement sur consigne non lisse, pas pour une mesure fine.
    """
    A = math.radians(amp_deg)
    T = 1.0 / freq_hz
    state = {"t": 0.0, "th": 0.0, "thd": 0.0, "u": 0.0}

    def ref(t):
        phi = (t / T) % 1.0
        if phi < 0.25:
            psi_r, psi_dot_r = A * 4 * phi, 4 * A / T
        elif phi < 0.75:
            psi_r, psi_dot_r = A * (2 - 4 * phi), -4 * A / T
        else:
            psi_r, psi_dot_r = A * (4 * phi - 4), 4 * A / T
        g = _fade(t)
        if not feedforward:
            return np.array([0.0, 0.0, g * psi_r, 0.0]), 0.0
        u_r = g * (-cfg.A21 * psi_r - cfg.A22 * psi_dot_r) / cfg.B2
        dt = max(t - state["t"], 0.0)
        state["th"] += state["thd"] * dt + 0.5 * state["u"] * dt ** 2
        state["thd"] += 0.5 * (state["u"] + u_r) * dt      # trapeze
        state["t"], state["u"] = t, u_r
        return (np.array([state["th"], state["thd"],
                          g * psi_r, g * psi_dot_r]), u_r)
    return ref


# =========================================================================== #
#  Verification d'admissibilite
# =========================================================================== #

def check_envelope(amp_deg, freq_hz):
    """Verifie que la trajectoire de reference tient dans les limites du
    banc. Le calcul est EXACT (pas de table interpolee) et porte sur la
    trajectoire elle-meme, pas sur une approximation."""
    u_max, thd_max, th_max = cfg.reference_effort(amp_deg, freq_hz)
    a_max = cfg.max_amplitude_deg(freq_hz)

    print("\n  effort demande au cerceau (trajectoire de reference seule) :")
    print("    max |u_r|         = {:7.1f} rad/s^2   (borne {:.0f})"
          .format(u_max, cfg.U_SLIP_MAX))
    print("    max |theta_dot_r| = {:7.2f} rad/s     (borne {:.0f})"
          .format(thd_max, cfg.THETA_DOT_MAX))
    print("    max |theta_r|     = {:7.2f} rad       ({:.2f} tour)"
          .format(th_max, th_max / (2 * math.pi)))
    print("    amplitude max a {:.2f} Hz : {:.1f} deg".format(freq_hz, a_max))

    ok = True
    if thd_max > cfg.THETA_DOT_MAX:
        print("\n[!!] {:.1f} deg a {:.2f} Hz demande {:.1f} rad/s au cerceau, "
              "soit {:.1f}x la limite.".format(amp_deg, freq_hz, thd_max,
                                               thd_max / cfg.THETA_DOT_MAX))
        print("     Le banc NE PEUT PAS produire ce mouvement : la balle ne")
        print("     suivra pas, quel que soit le regulateur.")
        ok = False
    elif thd_max > 0.7 * cfg.THETA_DOT_MAX:
        print("\n[!]  {:.0f} % de la vitesse disponible : peu de marge pour la"
              " correction d'ecart.".format(100 * thd_max / cfg.THETA_DOT_MAX))
    if u_max > cfg.U_SLIP_MAX:
        print("\n[!!] max |u_r| depasse la borne de non-glissement : la balle "
              "glisserait, hors modele.")
        ok = False
    return ok


def print_envelope_table():
    print("\n  enveloppe (theta_dot_max = {:.0f} rad/s, resonance a {:.3f} Hz)"
          .format(cfg.THETA_DOT_MAX, cfg.F_N_HZ))
    print("    {:>8s}  {:>12s}".format("f [Hz]", "A_max [deg]"))
    for f in (0.1, 0.2, 0.5, 0.8, 1.0, 1.2, cfg.F_N_HZ, 1.5, 2.0, 3.0):
        a = cfg.max_amplitude_deg(f)
        flag = "  <- hors domaine lineaire" if a > 30.0 else ""
        print("    {:8.3f}  {:12.1f}{}".format(f, a, flag))
    print("\n  Ces valeurs sont une borne SUPERIEURE issue du modele")
    print("  linearise autour de psi = 0. Au-dela de ~30 deg le modele")
    print("  petit-angle n'est plus valable et la borne n'a plus de sens :")
    print("  pres de la resonance, la contrainte reelle devient la")
    print("  non-linearite de sin(psi), pas la vitesse du cerceau.")
    print("\n  L'amplitude maximale passe par un maximum a la resonance :")
    print("  c'est le point de fonctionnement ou le cerceau a le moins de")
    print("  travail a fournir, donc le meilleur pour un essai de poursuite.")


# =========================================================================== #
#  Programme principal
# =========================================================================== #

def main():
    p = argparse.ArgumentParser(
        description="T2 : poursuite d'une reference variable")
    p.add_argument("--ref", choices=["sine", "step", "triangle"],
                   default=cfg.REF_TYPE)
    p.add_argument("--amp", type=float, default=cfg.REF_AMP_DEG,
                   help="amplitude de psi_ref [deg]")
    p.add_argument("--freq", type=float, default=cfg.REF_FREQ_HZ,
                   help="frequence de psi_ref [Hz]")
    p.add_argument("--gain", default=cfg.DEFAULT_GAIN,
                   choices=sorted(cfg.GAINS))
    p.add_argument("--estimator", choices=["derivative", "kalman", "ekf"],
                   default="derivative")
    p.add_argument("--no-feedforward", action="store_true",
                   help="desactive u_ff et x_ref complet : reproduit la v1, "
                        "pour comparaison chiffree")
    p.add_argument("--duration", type=float, default=20.0)
    p.add_argument("--trial", type=int, default=1)
    p.add_argument("--outdir", default=cfg.OUTDIR)
    p.add_argument("--envelope", action="store_true",
                   help="affiche l'enveloppe et quitte")
    p.add_argument("--force", action="store_true",
                   help="lance meme hors enveloppe")
    args = p.parse_args()

    if args.envelope:
        print_envelope_table()
        return

    if args.duration > cfg.MAX_DURATION_S:
        raise SystemExit("duree superieure au plafond de {:.0f} s"
                         .format(cfg.MAX_DURATION_S))

    bc.validate_config(cfg)

    K = bc.pick_gain(cfg, args.gain)
    ff = not args.no_feedforward

    print("=" * 70)
    print("  T2 v3 -- poursuite de reference, cerceau {}".format(cfg.HOOP))
    print("=" * 70)
    print("  reference    : {}, {:.1f} deg a {:.2f} Hz"
          .format(args.ref, args.amp, args.freq))
    print("  feedforward  : {}".format(
        "OUI (x_ref complet + u_ff)" if ff
        else "NON -- mode v1, x_ref = [0,0,psi_ref,0]"))
    print("  gains ({})  : K = [{:.4f}, {:.4f}, {:.4f}, {:.4f}]"
          .format(args.gain, *K))
    print("  theta        : {}".format(
        "encodeur" if cfg.USE_MEASURED_THETA else "integration de u"))
    print("  psi_dot      : {}".format(
        {"kalman": "Kalman lineaire + compensation de retard",
         "ekf": "Kalman etendu (sin psi) + compensation de retard"}
        .get(args.estimator,
             "derivee filtree {:.0f} Hz".format(cfg.PSIDOT_FILTER_HZ))))
    print("  ODrive       : consigne {:.0f} Hz, encodeur {:.0f} Hz, "
          "erreurs {:.0f} Hz".format(cfg.LOOP_HZ / cfg.COM_PERIOD,
                                     cfg.LOOP_HZ / cfg.ENC_PERIOD,
                                     cfg.LOOP_HZ / cfg.ERR_PERIOD))
    print("  modele       : f_n = {:.4f} Hz, zeta = {:.3f}, B2 = {:+.3f}"
          .format(cfg.F_N_HZ, cfg.ZETA, cfg.B2))

    A = np.array([[0, 1, 0, 0], [0, 0, 0, 0],
                  [0, 0, 0, 1], [0, 0, cfg.A21, cfg.A22]], float)
    B = np.array([0.0, 1.0, 0.0, cfg.B2])
    lam = np.linalg.eigvals(A - np.outer(B, K))
    print("  poles boucle fermee (modele) : max Re = {:+.3f}"
          .format(float(lam.real.max())))
    if lam.real.max() > 0:
        raise SystemExit(
            "\n[!!] boucle fermee instable dans le modele. Verifier le signe "
            "de B2 sur le banc avant de lancer T2 : avec le mauvais signe, "
            "u_ff envoie le cerceau a l'envers a pleine amplitude.")

    if args.ref == "sine":
        if not check_envelope(args.amp, args.freq) and not args.force:
            raise SystemExit("\nannule (utiliser --force pour passer outre).")
        reference = sine_reference(args.amp, args.freq, feedforward=ff)
    elif args.ref == "step":
        print("\n  [!] echelon : consigne NON ADMISSIBLE (voir en-tete). "
              "Erreur statique et derive de vitesse attendues.")
        reference = step_reference(args.amp, feedforward=ff)
    else:
        print("\n  [!] triangle : psi_ddot_r impulsionnelle, feedforward "
              "partiel et integration numerique. Essai qualitatif.")
        # Le triangle n'a pas d'enveloppe analytique, mais son harmonique
        # fondamentale a 8/pi^2 = 0.81 fois son amplitude : la sinusoide de
        # meme amplitude est une borne inferieure utile de l'effort.
        check_envelope(args.amp, args.freq)
        reference = triangle_reference(args.amp, args.freq, feedforward=ff)

    if cfg.REF_FADE_S > 0.0:
        print("\n  fondu de {:.1f} s : le feedforward n'y est exact qu'a "
              "l'ordre 0 en gpoint.".format(cfg.REF_FADE_S))
        print("  Les {:.0f} premieres secondes ne sont pas exploitables pour "
              "chiffrer la".format(cfg.REF_FADE_S))
        print("  poursuite ; analyse_tracking() n'utilise que la seconde "
              "moitie de l'essai.")
        if args.duration < 4.0 * cfg.REF_FADE_S:
            print("  [!] duree de {:.0f} s pour un fondu de {:.1f} s : la "
                  "seconde moitie contient".format(args.duration,
                                                   cfg.REF_FADE_S))
            print("      encore du transitoire. Allonger l'essai.")

    stopper = bc.Stopper()
    picam2, det, roi, centre = bc.setup_camera(cfg)

    estimator, n_delay = build_estimator(args.estimator, cfg)
    if estimator is not None:
        print("[est] compensation de retard sur {} pas ({:.0f} ms)"
              .format(n_delay, 1e3 * n_delay * cfg.DT_NOMINAL))
        print("[est] sigma_psi = {:.3f} deg, q_accel = {:.1f}"
              .format(cfg.SIGMA_PSI_DEG, cfg.Q_ACCEL))

    odrv0 = ax = None
    try:
        odrv0, ax = bc.connect_and_prepare(cfg)
        # La balle part TOUJOURS au repos au fond en T2 : la reference part
        # de zero et monte par le fondu.
        psi_offset, _, arm_stats = bc.arm(picam2, det, roi, centre, cfg, 0.0)

        log = bc.control_loop(
            ax, picam2, det, roi, centre, psi_offset,
            K, reference, args.duration, stopper, cfg,
            estimator=estimator, tag="t2")

        path = bc.save_and_summarise(
            log, cfg, args.outdir,
            "T2v3_{}_{}_A{:04.1f}deg_f{:05.2f}Hz_{}{}_{}_t{:02d}_{{stamp}}.npz"
            .format(cfg.HOOP, args.ref, args.amp, args.freq, args.gain,
                    "" if ff else "_noff", args.estimator, args.trial),
            {"experiment": "T2", "version": 3, "hoop": cfg.HOOP,
             "trial": args.trial, "ref_type": args.ref,
             "ref_amp_deg": args.amp, "ref_freq_hz": args.freq,
             "feedforward": bool(ff), "ref_fade_s": cfg.REF_FADE_S,
             "estimator": args.estimator,
             "estimator_n_delay": int(n_delay),
             "sigma_psi_deg": cfg.SIGMA_PSI_DEG, "q_accel": cfg.Q_ACCEL,
             "arm_offset_deg": math.degrees(arm_stats["offset"]),
             "arm_osc_deg": math.degrees(arm_stats["osc_amp"]),
             "arm_noise_deg": math.degrees(arm_stats["sigma"]),
             "psi_offset_rad": float(psi_offset),
             "roi": [int(v) for v in roi]},
            K, args.gain)

        if args.ref == "sine":
            analyse_tracking(log, args.freq)
        print("\n[out] {}".format(path))

    finally:
        bc.shutdown(ax, picam2, cfg)


# =========================================================================== #
#  Analyse de poursuite
# =========================================================================== #

def analyse_tracking(log, freq_hz):
    """Gain et dephasage par regression lineaire sur [sin, cos].

    Pour une reference sinusoidale, la moyenne de l'erreur ne dit rien : la
    reponse pertinente est le RAPPORT D'AMPLITUDE et le DEPHASAGE. On ajuste
    psi(t) ~ a sin(wt) + b cos(wt) par moindres carres -- probleme lineaire,
    donc solution exacte et sans reglage -- puis on compare a la reference
    ajustee de la meme facon. C'est l'analyse par onde sinusoidale
    classique en identification frequentielle [Ljung, "System
    Identification: Theory for the User", 2e ed., 1999].

    On n'utilise que les cycles apres etablissement (seconde moitie de
    l'essai), pour ne pas melanger transitoire et regime.
    """
    t = log["t"]
    n = t.size
    if n < 50:
        return
    sel = slice(int(n * 0.5), n)
    w = 2.0 * math.pi * freq_hz
    M = np.column_stack([np.sin(w * t[sel]), np.cos(w * t[sel]),
                         np.ones(t[sel].size)])

    def fit(y):
        c, *_ = np.linalg.lstsq(M, y[sel], rcond=None)
        return math.hypot(c[0], c[1]), math.atan2(c[1], c[0]), c[2]

    amp_y, ph_y, off_y = fit(log["psi"])
    amp_r, ph_r, _ = fit(log["psi_ref"])
    if amp_r < 1e-6:
        return
    err = log["psi"][sel] - log["psi_ref"][sel]

    print("\n  --- poursuite a {:.2f} Hz (seconde moitie de l'essai) ---"
          .format(freq_hz))
    print("    amplitude reference : {:6.2f} deg".format(math.degrees(amp_r)))
    print("    amplitude balle     : {:6.2f} deg".format(math.degrees(amp_y)))
    print("    rapport d'amplitude : {:6.3f}   ({:+.2f} dB)"
          .format(amp_y / amp_r, 20 * math.log10(amp_y / amp_r)))
    dphi = math.degrees(ph_y - ph_r)
    dphi = (dphi + 180.0) % 360.0 - 180.0
    print("    dephasage           : {:+6.1f} deg   (soit {:+.1f} ms)"
          .format(dphi, 1e3 * dphi / 360.0 / freq_hz))
    print("    biais residuel      : {:+6.2f} deg".format(math.degrees(off_y)))
    print("    erreur              : RMS {:.2f} deg, crete {:.2f} deg"
          .format(math.degrees(float(np.sqrt(np.mean(err ** 2)))),
                  math.degrees(float(np.max(np.abs(err))))))
    print("    a comparer a sigma_psi = {:.2f} deg (bruit de mesure)"
          .format(cfg.SIGMA_PSI_DEG))


if __name__ == "__main__":
    main()
