#!/usr/bin/env python3
"""
t4_run_v3.py -- T4 : pompage, largage, vol libre, rattrapage a l'interieur.

    python3 t4_run_v3.py --dry-run
    python3 t4_run_v3.py --trial 1
    python3 t4_run_v3.py --plan plans/T4d_flying_ball_plan.mat

AVERTISSEMENT -- LE PLAN COURANT N'EST PAS VOLABLE (28/08)
==========================================================
Ce fichier est le pilote. Il est termine et verifie, mais le plan qu'il
lit ne l'est pas : avec constrain_theta_f = true la simulation MATLAB sort
la sequence de modes "rolling_out -> free_fall -> rolling_out", c'est-a-
dire que la balle retombe sur le cerceau EXTERIEUR au lieu d'etre
rattrapee a l'interieur. Avec constrain_theta_f = false elle est rattrapee
en simulation mais passe dans la matiere du cerceau interieur sur le banc.
Voir l'en-tete de scenarios/closed_loop/T4/T4d_flying_ball.m.

--dry-run REFUSE de valider un plan dont la sequence de modes a echoue,
quand l'export le dit. Il n'y a pas de raison d'apprendre cela avec une
bille lancee a 0.67 m/s.

LES DEUX MODIFICATIONS DE bench_common_v3.py QUE T4 IMPOSE
==========================================================
Les deux sont FAITES dans bench_common_v3.py, et elles sont de meme
nature que les deux que T3 avait imposees.

1. LE RAYON EST MESURE ET JOURNALISE. control_loop calculait
   psi = atan2(u_px - cx, v_px - cy) et jetait la norme. r_px =
   hypot(des deux memes composantes) ne coute pas une ligne de traitement
   d'image, et c'est lui qui porte le test de phase. Il est desormais
   journalise pour tous les essais, T1/T2/T3 compris.

2. LA LOI DE COMMANDE PEUT DEPENDRE DE L'ETAT. T3 avait besoin d'un gain
   fonction du TEMPS (gain_schedule). T4 a besoin d'une loi fonction de
   L'ETAT MESURE : le passage du vol au rattrapage se decide sur le rayon,
   pas sur l'horloge, parce que l'instant du contact depend du vol reel et
   non du plan. control_loop accepte donc un phase_hook(t, x, r_px) qui
   renvoie (x_ref, u_ff, K, phase), avec K = None pour la boucle ouverte.
   phase_hook=None laisse T1/T2/T3 inchanges au bit pres.

LES TROIS PHASES
================
  1  t <= Tf                    TVLQR le long du plan, comme T3.
  2  t > Tf, r > seuil          u = 0. Rien de ce que fait le moteur ne
                                change ou la balle atterrit. En T4d la
                                branche "speed hold" du controleur MATLAB
                                degenere exactement a cela (kp_hold =
                                kd_hold = 0), donc ce n'est pas une
                                simplification : c'est la meme loi.
  3  verrou des r <= seuil      LQR stationnaire K_catch, cible theta_park.

VERROU, ET PAS TEST INSTANTANE. En simulation le test est |r - R| < tol ET
|r_dot| < tol. Sur le banc r_dot serait derive d'un rayon bruite a 50 Hz :
inutilisable pour une decision binaire. Une fois la balle rattrapee elle ne
repart pas, et si elle repart l'essai est rate de toute facon.

LE REPLIAGE, QUI EST LE PIEGE DE LA PHASE 3
===========================================
UNWRAP_PSI est vrai (la balle monte a 120 deg puis vole), donc psi
accumule. La cible de la phase 3 est psi = 0 : regler sur un psi accumule
a -400 deg commanderait une correction absurde. theta accumule aussi, et
seul theta MODULO 2*pi est physique -- le trou est ou il est quel que soit
le nombre de tours faits pour y arriver. Les deux sont donc replies au
moment de former l'erreur, et SEULEMENT la : l'estimateur continue de
voir le signal deroule, sinon psi_dot sauterait au repliage.
"""

import argparse
import json
import math
import os

import numpy as np
from scipy.io import loadmat

from common import bench_common as bc
from tasks import t4_config as cfg
from common import hoop_geometry as geo


TWO_PI = 2.0 * math.pi


def fold(angle, centre):
    """Ramene angle dans (centre - pi, centre + pi]."""
    return angle - TWO_PI * round((angle - centre) / TWO_PI)


# ----------------------------------------------------------------------
# Plan
# ----------------------------------------------------------------------
def load_plan(path):
    """Lit le .mat ecrit par export_t4_for_bench.m.

    La conversion de convention de signe (psi_banc = -psi_MATLAB, sur psi,
    psi_dot, K3, K4 ET K_catch(3), K_catch(4)) est faite DANS L'EXPORT.
    Un seul point de conversion, sinon elle finit par etre appliquee deux
    fois ou zero fois.
    """
    if not os.path.exists(path):
        raise SystemExit(
            "plan introuvable : {}\nGenerer d'abord cote MATLAB :\n"
            "  export_t4_for_bench(@T4d_flying_ball, "
            "'<...>/Python/t1_t2/plans/T4d_flying_ball_plan.mat')".format(path))

    m = loadmat(path, squeeze_me=True)
    plan = {
        "t":          np.asarray(m["t"], float).ravel(),
        "x_ref":      np.asarray(m["x_ref"], float),
        "u_ff":       np.asarray(m["u_ff"], float).ravel(),
        "K_traj":     np.asarray(m["K_traj"], float),
        "Tf":         float(m["Tf"]),
        "K_catch":    np.asarray(m["K_catch"], float).ravel(),
        "R_out":      float(m["R_out"]),
        "R_in":       float(m["R_in"]),
        "theta_park": float(m["theta_park"]),
    }
    n = plan["t"].size
    if plan["x_ref"].shape != (n, 4) or plan["K_traj"].shape != (n, 4):
        raise SystemExit(
            "formes incoherentes dans le plan : t {}, x_ref {}, K_traj {}"
            .format(plan["t"].shape, plan["x_ref"].shape, plan["K_traj"].shape))
    if plan["K_catch"].size != 4:
        raise SystemExit("K_catch doit avoir 4 termes, en a {}"
                         .format(plan["K_catch"].size))
    return plan


def scale_px_per_m(plan):
    """Echelle image, deduite du rayon d'orbite EXTERIEUR.

    R_out metres valent HOOP_RADIUS_PX pixels. On ne prend pas
    cfg.R_EFF_M : le plan porte sa propre valeur de R_out, et si les deux
    divergent c'est la geometrie MATLAB qui fait foi puisque c'est elle
    qui a produit R_in.
    """
    return cfg.HOOP_RADIUS_PX / plan["R_out"]


def check_plan_against_config(plan):
    """Refuse un plan que la config du banc ecreterait, ou que la camera
    ne suivrait pas. Un plan tronque par une saturation n'est plus le
    plan : le TVLQR corrigerait une trajectoire qu'il ne suit pas."""
    problems = []
    px_per_m = scale_px_per_m(plan)

    thd = float(np.max(np.abs(plan["x_ref"][:, 1])))
    if thd > cfg.THETA_DOT_MAX:
        problems.append(
            "le plan demande |theta_dot| = {:.2f} rad/s, THETA_DOT_MAX = {:.2f}"
            .format(thd, cfg.THETA_DOT_MAX))

    umax = float(np.max(np.abs(plan["u_ff"])))
    if umax > cfg.U_SLIP_MAX:
        problems.append(
            "le plan demande |u| = {:.2f} rad/s^2, U_SLIP_MAX = {:.2f}"
            .format(umax, cfg.U_SLIP_MAX))

    u_hi = (cfg.TAU_MAX - cfg.B_MOTEUR * thd) / cfg.I_TOTAL
    if umax > u_hi:
        problems.append(
            "le plan demande |u| = {:.2f} rad/s^2, couple disponible a "
            "{:.1f} rad/s = {:.2f}".format(umax, thd, u_hi))

    psi_span = float(np.max(np.abs(plan["x_ref"][:, 2])))
    if math.degrees(psi_span) > cfg.ABORT_PSI_DEG:
        problems.append(
            "le plan atteint psi = {:.0f} deg, ABORT_PSI_DEG = {:.0f}"
            .format(math.degrees(psi_span), cfg.ABORT_PSI_DEG))

    # --- Fenetre de suivi : roulage ET vol ------------------------------
    psid_max = float(np.max(np.abs(plan["x_ref"][:, 3])))
    d_roll = cfg.HOOP_RADIUS_PX * psid_max / cfg.LOOP_HZ
    # Vol : composante tangentielle au largage, plus la chute radiale
    # moyenne sur la traversee R_out -> R_in.
    v_tan = plan["R_out"] * abs(float(plan["x_ref"][-1, 3]))
    d_fly = (v_tan / cfg.LOOP_HZ) * px_per_m
    d_px = max(d_roll, d_fly)
    if d_px > 0.5 * cfg.TRACKING_WIN:
        problems.append(
            "la balle se deplace de {:.0f} px entre deux images ({:.0f} en "
            "roulage, {:.0f} en vol),\n       au-dela de la demi-fenetre "
            "TRACKING_WIN/2 = {:.0f} px : le detecteur la perdra. Porter\n"
            "       TRACKING_WIN a {:.0f} au moins."
            .format(d_px, d_roll, d_fly, 0.5 * cfg.TRACKING_WIN,
                    4 * math.ceil(d_px / 2)))

    # --- Le seuil de rattrapage doit etre separable du bruit -------------
    r_in_px = plan["R_in"] * px_per_m
    catch_px = (plan["R_in"] + cfg.CATCH_TOL_M) * px_per_m
    if catch_px >= 0.9 * cfg.HOOP_RADIUS_PX:
        problems.append(
            "le seuil de rattrapage ({:.0f} px) est trop proche de l'orbite "
            "exterieure ({:.0f} px) :\n       la phase 3 se declencherait "
            "avant le largage".format(catch_px, cfg.HOOP_RADIUS_PX))
    if r_in_px < 20.0:
        problems.append(
            "R_in ne fait que {:.0f} px : le rayon mesure y est domine par "
            "le bruit de\n       detection, le test de phase serait un tirage "
            "au sort. Verifier HOOP_RADIUS_PX.".format(r_in_px))

    # --- Stabilite du gain de rattrapage --------------------------------
    # Detecte une erreur de SIGNE dans l'export, qui est le mode de defaut
    # redoute (psi_banc = -psi_MATLAB). Le test a un sens ici parce que la
    # phase 3 EST une regulation stationnaire au fond, la ou le modele
    # linearise est exactement valable. ATTENTION : le modele du banc
    # charge par t1_config depend de T1_HOOP, et la phase 3 se passe sur la
    # face INTERNE du cerceau interieur -- voir la note plus bas.
    A = np.array([[0, 1, 0, 0], [0, 0, 0, 0],
                  [0, 0, 0, 1], [0, 0, cfg.A21, cfg.A22]], float)
    B = np.array([0.0, 1.0, 0.0, cfg.B2])
    lam = np.linalg.eigvals(A - np.outer(B, plan["K_catch"]))
    if lam.real.max() > 0:
        problems.append(
            "K_catch donne une boucle fermee INSTABLE dans le modele "
            "{} (max Re = {:+.3f}).\n       Verifier la convention de signe "
            "de l'export : la symetrie psi_banc = -psi_MATLAB\n       inverse "
            "le QUADRUPLET complet ET u_ff, auquel cas K est INCHANGE."
            .format(cfg.HOOP, float(lam.real.max())))

    # --- Le plan tient-il ce qu'il promet en simulation ? ---------------
    thd_end = float(plan["x_ref"][-1, 1])
    if abs(thd_end) < 0.5:
        problems.append(
            "le plan finit a theta_dot = {:+.2f} rad/s : le cerceau serait "
            "a l'arret pendant le vol.\n       C'est la rotation du cerceau "
            "qui dissipe l'energie de la balle a l'impact ; a l'arret elle "
            "est\n       renvoyee par le trou.".format(thd_end))

    return problems


# ----------------------------------------------------------------------
# Machine a trois phases
# ----------------------------------------------------------------------
def build_phase_hook(plan, verbose=True):
    """Retourne (hook, state). hook(t, x, r_px) -> (x_ref, u_ff, K, phase).

    state est un dict inspecte apres l'essai : instant et rayon du
    basculement, phase atteinte. Il n'entre dans aucune decision.
    """
    t_grid, X, U, Kt = plan["t"], plan["x_ref"], plan["u_ff"], plan["K_traj"]
    Tf = plan["Tf"]
    K_catch = plan["K_catch"]
    theta_park = plan["theta_park"]
    catch_px = (plan["R_in"] + cfg.CATCH_TOL_M) * scale_px_per_m(plan)

    state = {"phase": 1, "t_release": None, "t_catch": None,
             "r_catch_px": None, "catch_px": catch_px}

    def hook(t, x, r_px):
        # --- verrou : une fois en phase 3, on y reste -------------------
        if state["phase"] == 3:
            return _catch_ref(x), 0.0, K_catch, 3

        if t <= Tf:
            x_ref = np.array([np.interp(t, t_grid, X[:, j]) for j in range(4)])
            K = np.array([np.interp(t, t_grid, Kt[:, j]) for j in range(4)])
            return x_ref, float(np.interp(t, t_grid, U)), K, 1

        if state["phase"] == 1:
            state["phase"] = 2
            state["t_release"] = t
            if verbose:
                print("\n[t4] phase 2 -- vol libre, u = 0 (t = {:.3f} s, "
                      "r = {:.0f} px)".format(t, r_px))

        if r_px <= catch_px:
            state["phase"] = 3
            state["t_catch"] = t
            state["r_catch_px"] = r_px
            if verbose:
                print("\n[t4] phase 3 -- rattrapage verrouille (t = {:.3f} s, "
                      "r = {:.0f} px <= {:.0f})".format(t, r_px, catch_px))
            return _catch_ref(x), 0.0, K_catch, 3

        # Vol libre. x_ref est celui du dernier noeud du plan : il ne sert
        # qu'au journal, puisque K = None coupe le retour d'etat.
        return X[-1], 0.0, None, 2

    def _catch_ref(x):
        """Cible de la phase 3, avec les DEUX repliages.

        theta est replie sur l'angle de garage : seul theta modulo 2*pi est
        physique. psi est replie sur 0 : il a accumule pendant le pompage
        et le vol, et la cible du rattrapage est le fond du cerceau
        interieur, pas "zero tours plus tard".
        """
        # x[0] - fold(x[0], c) vaut 2*pi * (nombre de tours accumules) :
        # la cible est donc l'angle de garage decale du meme nombre de
        # tours, ce qui rend l'erreur x - x_ref egale a l'erreur repliee.
        return np.array([theta_park + (x[0] - fold(x[0], theta_park)),
                         0.0,
                         x[2] - fold(x[2], 0.0),
                         0.0])

    return hook, state


# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="T4 : pompage, largage, vol libre, rattrapage")
    p.add_argument("--plan", default=cfg.PLAN_FILE)
    p.add_argument("--estimator", choices=["derivative", "ekf"],
                   default=getattr(cfg, "DEFAULT_ESTIMATOR", "derivative"),
                   help="derivee filtree (defaut) : l'EKF porte le modele "
                        "du pendule exterieur, faux pendant le vol et sur "
                        "la face interne")
    p.add_argument("--trial", type=int, default=1)
    p.add_argument("--outdir", default=cfg.OUTDIR)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    plan = load_plan(args.plan)
    px_per_m = scale_px_per_m(plan)
    duration = plan["Tf"] + cfg.HOLD_AFTER_CATCH_S

    print("=" * 70)
    print("  T4 v3 -- balle volante, cerceau {}".format(cfg.HOOP))
    print("=" * 70)
    print("  plan       : {}".format(os.path.basename(args.plan)))
    print("  Tf         : {:.4f} s  (+ {:.1f} s de maintien)"
          .format(plan["Tf"], cfg.HOLD_AFTER_CATCH_S))
    print("  noeuds     : {} ({:.1f} Hz moyen, boucle a {:.0f} Hz)"
          .format(plan["t"].size, (plan["t"].size - 1) / plan["Tf"], cfg.LOOP_HZ))
    print("  largage    : psi = {:+.1f} deg, psi_dot = {:+.2f} rad/s, "
          "theta_dot = {:+.2f} rad/s"
          .format(math.degrees(plan["x_ref"][-1, 2]), plan["x_ref"][-1, 3],
                  plan["x_ref"][-1, 1]))
    print("  theta final: {:+.2f} deg mod 360   (place le trou)"
          .format(math.degrees(plan["x_ref"][-1, 0]) % 360.0))
    print("  max |u_ff| : {:.2f} rad/s^2   (U_SLIP_MAX = {:.0f})"
          .format(np.max(np.abs(plan["u_ff"])), cfg.U_SLIP_MAX))
    print("  max |th_d| : {:.2f} rad/s     (THETA_DOT_MAX = {:.0f})"
          .format(np.max(np.abs(plan["x_ref"][:, 1])), cfg.THETA_DOT_MAX))
    print("  max |psi_d|: {:.2f} rad/s"
          .format(np.max(np.abs(plan["x_ref"][:, 3]))))
    print("  echelle    : {:.1f} px/m  ({:.1f} px pour R_out = {:.1f} mm)"
          .format(px_per_m, cfg.HOOP_RADIUS_PX, 1e3 * plan["R_out"]))
    print("  rayons     : R_out {:.0f} px, R_in {:.0f} px, seuil de "
          "rattrapage {:.0f} px"
          .format(plan["R_out"] * px_per_m, plan["R_in"] * px_per_m,
                  (plan["R_in"] + cfg.CATCH_TOL_M) * px_per_m))
    print("  theta_park : {:+.4f} rad".format(plan["theta_park"]))
    print("  K_catch    : [{:+.4f}, {:+.4f}, {:+.4f}, {:+.4f}]"
          .format(*plan["K_catch"]))
    print("  psi_dot    : {}".format(
        "Kalman etendu" if args.estimator == "ekf"
        else "derivee filtree {:.0f} Hz".format(cfg.PSIDOT_FILTER_HZ)))

    problems = check_plan_against_config(plan)
    problems += geo.check_geometry(cfg)
    if problems:
        print("\n[!!] Le plan et la configuration du banc sont incompatibles :")
        for q in problems:
            print("     - {}".format(q))
        raise SystemExit("corriger t4_config_v3.py ou replanifier cote MATLAB")
    print("\n  plan compatible avec la configuration du banc")

    print("\n  [!] CE CONTROLE NE DIT PAS QUE LE PLAN RATTRAPE LA BALLE.")
    print("      Il verifie que le banc peut l'executer, pas que la balle")
    print("      passe le trou ni qu'elle se pose a l'interieur. Ces deux")
    print("      questions se tranchent cote MATLAB : sequence de modes")
    print("      attendue dans le run, et t4d_gap_clearance.m pour le trou.")
    print("      Au 28/08 le plan T4d echoue la premiere. NE PAS LANCER.")

    if args.dry_run:
        print("\n  --dry-run : rien n'a ete envoye au materiel")
        return

    bc.validate_config(cfg)

    stopper = bc.Stopper()
    picam2, det, roi, centre = bc.setup_camera(cfg)

    if args.estimator == "ekf":
        from common import ekf as est_mod
        estimator, n_delay = est_mod.build_estimator(cfg, cfg.DT_NOMINAL)
    else:
        estimator, n_delay = None, 0

    hook, hook_state = build_phase_hook(plan)

    ax = None
    try:
        odrv0, ax = bc.connect_and_prepare(cfg)
        psi_offset, psi_init_meas, arm_stats = bc.arm(
            picam2, det, roi, centre, cfg, 0.0)

        log = bc.control_loop(
            ax, picam2, det, roi, centre, psi_offset,
            plan["K_traj"][0], None, duration, stopper, cfg,
            estimator=estimator, tag="t4", phase_hook=hook)

        print("\n  phase atteinte      : {}".format(hook_state["phase"]))
        print("  largage             : {}".format(
            "t = {:.3f} s".format(hook_state["t_release"])
            if hook_state["t_release"] else "jamais atteint"))
        print("  rattrapage          : {}".format(
            "t = {:.3f} s a r = {:.0f} px".format(hook_state["t_catch"],
                                                  hook_state["r_catch_px"])
            if hook_state["t_catch"] else "NON DETECTE"))

        bc.save_and_summarise(
            log, cfg, args.outdir,
            "T4v3_{}_fly_{}_t{:02d}_{{stamp}}.npz".format(
                cfg.HOOP, args.estimator, args.trial),
            {"experiment": "T4", "version": 3, "hoop": cfg.HOOP,
             "trial": args.trial,
             "plan_file": os.path.basename(args.plan),
             "Tf": plan["Tf"],
             "psi_release_deg": math.degrees(plan["x_ref"][-1, 2]),
             "theta_dot_release": float(plan["x_ref"][-1, 1]),
             "theta_final_deg_mod360": math.degrees(plan["x_ref"][-1, 0]) % 360.0,
             "K_catch": [float(v) for v in plan["K_catch"]],
             "theta_park": plan["theta_park"],
             "R_out_m": plan["R_out"], "R_in_m": plan["R_in"],
             "catch_tol_m": cfg.CATCH_TOL_M,
             "catch_threshold_px": hook_state["catch_px"],
             "px_per_m": px_per_m,
             "phase_reached": hook_state["phase"],
             "t_release_s": hook_state["t_release"],
             "t_catch_s": hook_state["t_catch"],
             "r_catch_px": hook_state["r_catch_px"],
             "arm_offset_deg": math.degrees(arm_stats["offset"]),
             "arm_noise_deg": math.degrees(arm_stats["sigma"]),
             "ref_type": "t4_three_phase", "estimator": args.estimator,
             "estimator_n_delay": int(n_delay),
             "tracking_win": int(cfg.TRACKING_WIN),
             "hoop_radius_px": cfg.HOOP_RADIUS_PX,
             "unwrap_psi": bool(cfg.UNWRAP_PSI),
             "psi_offset_rad": float(psi_offset),
             "roi": [int(v) for v in roi]},
            plan["K_traj"][0], "t4_three_phase")

    finally:
        bc.shutdown(ax, picam2, cfg)


if __name__ == "__main__":
    main()
