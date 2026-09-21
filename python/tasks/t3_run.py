#!/usr/bin/env python3
"""
t3_run_v3.py -- T3 : looping complet, trajectoire planifiee + TVLQR.

Meme architecture que t1_run_v3.py et t2_run_v3.py : toute la mecanique est
dans bench_common_v3.py, ce fichier ne fournit que le generateur de
reference -- ici lu dans un plan exporte par MATLAB plutot que calcule en
ligne.

    x_ref(t), u_ff(t), K(t)   <-  export_tvlqr_for_bench.m

LES DEUX MODIFICATIONS DE bench_common_v3.py QUE T3 IMPOSE
==========================================================
Les deux sont FAITES dans bench_common_v3.py.

1. GAIN VARIANT. control_loop() applique
        K_t = K if gain_schedule is None else gain_schedule(t_rel)
        u_raw = u_ff - gain_scale * float(K_t @ (x - x_ref))
   T1 et T2 passent gain_schedule=None et sont donc bit-identiques. Rien
   dans la chaine d'acquisition (camera, cadences, ODriveLink, journal)
   n'est touche : le patch est dans la loi de commande, pas le pipeline.

2. DEROULAGE DE psi. psi vient de atan2, donc replie dans (-pi, pi]. Sur
   un tour complet il saute de +180 a -180 deg. Ce saut entrerait dans le
   derivateur de psi_dot (bond de 2*pi en 20 ms = 314 rad/s de vitesse
   fantome) et dans l'erreur de suivi. L'accumulateur de control_loop()
   ramene chaque increment dans (-pi, pi] avant de le cumuler. Il est
   valide tant que la balle tourne de moins d'un demi-tour ENTRE DEUX
   DETECTIONS REUSSIES : c'est donc ABORT_LOST_FRAMES qui en fixe la
   marge, pas la cadence de boucle (0.16 tour pour 5 images perdues a
   1.6 rev/s, contre 0.5 admissible).

ESTIMATION DE PSI_DOT SUR UN TOUR COMPLET
=========================================
--estimator kalman reste INTERDIT. Le filtre d'estimator_v3 porte le
modele linearise autour de psi = 0 : le rappel y vaut A21*psi au lieu de
A21*sin(psi), et les deux sont de SIGNE OPPOSE au-dela de 180 deg,
c'est-a-dire au sommet du looping. Comme control_loop() remplace psi par
l'estimation, un filtre hors domaine ne degrade pas seulement psi_dot : il
corrompt la mesure. Le refus est porte par estimator_v3.build_estimator
(via UNWRAP_PSI) et double ici par la liste de choix.

--estimator ekf leve cette limite. ekf_v3.ExtendedKF porte le rappel
sinusoidal et relinearise la jacobienne (A21*cos(psi)) a chaque pas, donc
il reste valable sur un tour complet. En psi = 0 il coincide exactement
avec le filtre lineaire, ce qui rend la comparaison T1/T2 sans autre
variable.

    ekf         Kalman etendu + compensation de retard. CHOIX PAR DEFAUT
                depuis que sa validation sur T1 est acquise (prerequis pose
                par la version precedente de cet en-tete, et leve).
    derivative  derivee arriere + passe-bas, sans modele. Conserve pour
                reproduire les campagnes anterieures au changement, et comme
                repli si l'EKF devait etre mis en cause.

POURQUOI L'EKF PLUTOT QUE LA DERIVEE FILTREE, EN T3 PARTICULIEREMENT
====================================================================
build_estimator() renvoie n_delay = 0 pour "derivative" : la latence
camera n'est alors PAS compensee du tout, alors qu'elle l'etait sur un pas
en T1/T2. C'est une regression propre a T3, et elle porte sur la voie qui
compte -- mesure faite sur les huit essais du 10/08 (Analyse/
t3_delay_budget.py) : le filtre lui-meme n'ajoute que 10.6 ms mediane
(attendu 23.3), la cadence est saine (0 % de pas > 2.5 x median) et le
variateur suit (0.0-0.9 ms de retard, vel_enc = 0.95 * vel_cmd). Le retard
restant est donc essentiellement la latence camera non compensee.

n_delay = int(round(CAM_LATENCY_S / DT_NOMINAL)) = 1 pas. A 20 ms de pas,
TOUTE latence entre ~11 et ~29 ms donne le meme n_delay = 1 : l'incertitude
sur CAM_LATENCY_S (19.3 ms mesures a 18000 us de pose, T3 tournant a 9000)
ne change donc pas la compensation. La re-mesure reste souhaitable pour
citer un chiffre, elle n'est pas bloquante pour l'usage.

Usage
-----
    python3 t3_run_v3.py --dry-run          # verifie le plan, ne bouge pas
    python3 t3_run_v3.py --trial 1
    python3 t3_run_v3.py --estimator ekf --trial 1
    python3 t3_run_v3.py --plan plans/autre_plan.mat

Arret : Ctrl-C. La rampe de deceleration est appliquee dans tous les cas.
"""

import argparse
import math
import os

import numpy as np
from scipy.io import loadmat

from common import bench_common as bc
from tasks import t3_config as cfg
from common import hoop_geometry as geo

def build_estimator(name, cfg):
    """Retourne (estimateur, n_delay) ou (None, 0) pour la derivee filtree.

    'kalman' n'est pas propose : voir l'en-tete. Si le nom arrivait quand
    meme ici, estimator_v3.build_estimator le refuserait de toute facon sur
    UNWRAP_PSI -- la double garde est volontaire.
    """
    if name == "derivative":
        return None, 0
    if name == "ekf":
        from common import ekf as est_mod
        return est_mod.build_estimator(cfg, cfg.DT_NOMINAL)
    raise SystemExit("estimateur inconnu ou interdit en T3 : {}".format(name))


# ----------------------------------------------------------------------
# Chargement du plan
# ----------------------------------------------------------------------
def load_plan(path):
    """Lit le .mat ecrit par export_tvlqr_for_bench.m.

    Retour : dict avec t (K,), x_ref (K,4), u_ff (K,), K_traj (K,4), Tf.

    La conversion de convention de signe (psi_banc = -psi_MATLAB) est faite
    DANS L'EXPORT, pas ici : un seul point de conversion, sinon elle finit
    par etre appliquee deux fois ou zero fois.
    """
    if not os.path.exists(path):
        raise SystemExit(
            "plan introuvable : {}\n"
            "Generer d'abord cote MATLAB :\n"
            "  export_tvlqr_for_bench(@T3c_loop_no_slip, "
            "'<...>/Python/t1_t2/plans/T3c_loop_no_slip_plan.mat')".format(path))

    m = loadmat(path, squeeze_me=True)
    plan = {
        "t":      np.asarray(m["t"], float).ravel(),
        "x_ref":  np.asarray(m["x_ref"], float),
        "u_ff":   np.asarray(m["u_ff"], float).ravel(),
        "K_traj": np.asarray(m["K_traj"], float),
        "Tf":     float(m["Tf"]),
    }
    n = plan["t"].size
    if plan["x_ref"].shape != (n, 4) or plan["K_traj"].shape != (n, 4):
        raise SystemExit(
            "formes incoherentes dans le plan : t {}, x_ref {}, K_traj {}"
            .format(plan["t"].shape, plan["x_ref"].shape, plan["K_traj"].shape))
    return plan


def check_plan_against_config(plan):
    """Refuse un plan que la config du banc ecreterait.

    Un plan tronque par une saturation n'est plus le plan : le TVLQR
    corrigerait une trajectoire qu'il ne suit pas. Mieux vaut refuser de
    partir que produire un essai ininterpretable.
    """
    problems = []

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

    if not getattr(cfg, "UNWRAP_PSI", False) and math.degrees(psi_span) > 170.0:
        problems.append(
            "le plan depasse 170 deg mais UNWRAP_PSI est faux : psi se "
            "replierait et psi_dot exploserait au passage")

    # --- Fenetre de suivi du detecteur -----------------------------------
    # Le file pendant la pose n'est PAS la contrainte dominante : c'est le
    # deplacement ENTRE DEUX IMAGES, qui doit rester dans la demi-fenetre
    # de suivi. Sans HOOP_RADIUS_PX on ne peut pas trancher, et partir sans
    # avoir tranche revient a decouvrir le probleme au sommet du looping.
    psid_max = float(np.max(np.abs(plan["x_ref"][:, 3])))
    r_px = getattr(cfg, "HOOP_RADIUS_PX", None)
    if r_px is None:
        problems.append(
            "HOOP_RADIUS_PX n'est pas renseigne dans t3_config_v3.py : "
            "impossible de verifier\n       que le deplacement de la balle "
            "entre deux images ({:.1f} mm au pic) tient\n       dans la "
            "demi-fenetre de suivi TRACKING_WIN/2 = {:.0f} px. Mesurer le "
            "rayon\n       du cercle decrit par le centre de la balle, en "
            "pixels, et le figer."
            .format(1e3 * cfg.R_EFF_M * psid_max / cfg.LOOP_HZ,
                    cfg.TRACKING_WIN / 2.0))
    else:
        d_px = r_px * psid_max / cfg.LOOP_HZ
        if d_px > 0.5 * cfg.TRACKING_WIN:
            problems.append(
                "la balle se deplace de {:.0f} px entre deux images au pic, "
                "au-dela de la demi-fenetre\n       de suivi "
                "TRACKING_WIN/2 = {:.0f} px : le detecteur la perdra au "
                "sommet. Porter\n       TRACKING_WIN a {:.0f} au moins."
                .format(d_px, 0.5 * cfg.TRACKING_WIN, 4 * math.ceil(d_px / 2)))

    # --- La ROI doit contenir la BALLE ENTIERE ----------------------------
    # ROI_RADIUS_PX dimensionne la ROI ; HOOP_RADIUS_PX decrit l'orbite du
    # CENTRE de la balle. Les confondre donne une ROI tangente a l'orbite,
    # donc une balle rognee sur tout son parcours -- c'est ce qui a fait
    # echouer la campagne du 10/08.
    r_roi = getattr(cfg, "ROI_RADIUS_PX", getattr(cfg, "HOOP_RADIUS_PX", None))
    if r_px is not None and r_roi is not None:
        r_ball_px = 0.0125 / (cfg.R_EFF_M / r_px)     # 12.5 mm en pixels
        besoin = r_px + r_ball_px
        dispo = r_roi * (1.0 + cfg.ROI_MARGIN)
        if dispo < besoin:
            problems.append(
                "la ROI couvre {:.0f} px de rayon alors que la balle entiere "
                "en demande {:.0f}\n       ({:.0f} d'orbite + {:.0f} de rayon "
                "de balle) : elle sera rognee. Porter\n       ROI_RADIUS_PX a "
                "{:.0f} au moins.".format(dispo, besoin, r_px, r_ball_px,
                                          besoin / (1.0 + cfg.ROI_MARGIN)))

    # --- Etat tenu apres Tf ----------------------------------------------
    # Au-dela de Tf la reference gele x_ref[-1]. Si le plan y laisse une
    # vitesse de cerceau non nulle, le retour d'etat la maintient pendant
    # tout HOLD_AFTER_TF_S : le cerceau tourne a vitesse constante apres la
    # manoeuvre, ce qui n'est pas une regulation stationnaire.
    thd_end = float(plan["x_ref"][-1, 1])
    psid_end = float(plan["x_ref"][-1, 3])
    if abs(thd_end) > 0.5:
        problems.append(
            "le plan finit a theta_dot = {:+.2f} rad/s : la reference gelee "
            "apres Tf commanderait\n       cette vitesse pendant les {:.1f} s "
            "de maintien. Replanifier avec une fin au repos."
            .format(thd_end, cfg.HOLD_AFTER_TF_S))
    if abs(psid_end) > 1.0:
        problems.append(
            "le plan finit a psi_dot = {:+.2f} rad/s : la cible de maintien "
            "n'est pas un point\n       d'equilibre, le TVLQR regulerait vers "
            "un etat que la balle ne peut pas tenir."
            .format(psid_end))

    # --- Stabilite du gain aux DEUX EXTREMITES ---------------------------
    # Le plan part et arrive au fond du cerceau (psi = 0 et psi = 2*pi), donc
    # le modele linearise autour du fond y est exactement valable et le test
    # sur les valeurs propres de (A - B K) y a un sens rigoureux. Aux
    # instants intermediaires il n'en aurait PAS : les valeurs propres gelees
    # d'un systeme variant dans le temps ne caracterisent pas sa stabilite
    # [Khalil, "Nonlinear Systems", 3e ed., 2002, sect. 4.5]. Le test n'est
    # donc pas fait ailleurs, et il ne pretend pas valider le TVLQR : il
    # detecte une erreur de SIGNE dans l'export, qui est le mode de defaut
    # redoute (psi_banc = -psi_MATLAB). t1_run_v3 et t2_run_v3 font deja ce
    # controle avant tout mouvement ; T3 ne le faisait pas.
    A = np.array([[0, 1, 0, 0], [0, 0, 0, 0],
                  [0, 0, 0, 1], [0, 0, cfg.A21, cfg.A22]], float)
    B = np.array([0.0, 1.0, 0.0, cfg.B2])
    for label, Kx in (("initial K_traj[0]", plan["K_traj"][0]),
                      ("de maintien K_traj[-1]", plan["K_traj"][-1])):
        lam = np.linalg.eigvals(A - np.outer(B, Kx))
        if lam.real.max() > 0:
            problems.append(
                "le gain {} donne une boucle fermee INSTABLE dans le modele "
                "au fond\n       (max Re = {:+.3f}). Verifier la convention "
                "de signe de l'export : la\n       symetrie psi_banc = "
                "-psi_MATLAB inverse le QUADRUPLET complet\n       (theta, "
                "theta_dot, psi, psi_dot) ET u_ff, auquel cas K est "
                "INCHANGE.\n       N'inverser que psi et psi_dot en "
                "retournant K3 et K4 casse la coherence\n       avec "
                "B2 = {:+.3f}.".format(label, float(lam.real.max()), cfg.B2))

    return problems


def check_config_complete(cfg):
    """Verifie que t3_config expose tout ce que le pipeline consomme.

    T3 herite de t1_config. Une liste d'imports explicite se desynchronise
    en silence, et le trou ne se voit qu'au premier essai reel puisque
    --dry-run ne touche ni la camera ni l'ODrive. Ce controle est fait
    AVANT --dry-run, pour que le mode qui sert justement a se rassurer
    verifie aussi cela.
    """
    needed = (
        "THRESHOLD", "DOWNSAMPLE", "TRACKING_WIN", "BALL_SIZE", "ROI_MARGIN",
        "HOOP_CENTRE_PX", "USE_HOUGH_CENTRE", "SIZE", "FPS", "EXPOSURE_US",
        "GAIN_CAM", "CAM_LATENCY_S", "RAMP_DOWN_S",
        "ARM_SETTLE_S", "ARM_HOLD_S", "ARM_PSI_MAX_DEG", "ARM_OSC_MAX_DEG",
        "ARM_NOISE_MAX_DEG", "ARM_INIT_TOL_DEG", "F_N_HZ",
        "LOOP_HZ", "DT_NOMINAL", "DT_MAX",
        "COM_PERIOD", "ENC_PERIOD", "IQ_PERIOD", "ERR_PERIOD",
        "USE_MEASURED_THETA", "PSIDOT_FILTER_HZ", "GAIN_START", "GAIN_RAMP_S",
        "TAU_MAX", "B_MOTEUR", "I_TOTAL", "KT",
        "A21", "A22", "B2", "U_SLIP_MAX", "THETA_DOT_MAX",
        "ABORT_PSI_DEG", "ABORT_LOST_FRAMES", "ABORT_OVERRUN",
        "MAX_DURATION_S", "R_EFF_M", "SIGMA_PSI_DEG", "Q_ACCEL", 
        "TAU_VEL_S", "HOOP_RADIUS_PX"
    )
    missing = [n for n in needed if not hasattr(cfg, n)]
    if missing:
        raise SystemExit(
            "t3_config_v3.py est incomplet -- {} constante(s) manquante(s) :\n"
            "  {}\n"
            "Le pipeline planterait au premier acces (typiquement "
            "bc.setup_camera).\n"
            "Corriger par 'from t1_config_v3 import *' en tete de "
            "t3_config_v3.py."
            .format(len(missing), ", ".join(missing)))


def trajectory_reference(plan):
    """reference(t) -> (x_ref[4], u_ff) par interpolation lineaire.

    Au-dela de Tf on tient le dernier point : c'est exactement ce que fait
    tvlqr_controller.m cote simulation (branche t > Tf), donc les deux
    cotes regulent la meme chose apres la manoeuvre.
    """
    t, X, U, Tf = plan["t"], plan["x_ref"], plan["u_ff"], plan["Tf"]

    def ref(tau):
        if tau >= Tf:
            return X[-1].copy(), 0.0
        x = np.array([np.interp(tau, t, X[:, j]) for j in range(4)])
        return x, float(np.interp(tau, t, U))
    return ref


def gain_schedule(plan):
    """K(t) par interpolation lineaire, gele au dernier gain apres Tf."""
    t, Kt, Tf = plan["t"], plan["K_traj"], plan["Tf"]

    def sched(tau):
        if tau >= Tf:
            return Kt[-1]
        return np.array([np.interp(tau, t, Kt[:, j]) for j in range(4)])
    return sched


# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="T3 : looping complet, trajectoire + TVLQR")
    p.add_argument("--plan", default=cfg.PLAN_FILE,
                   help="fichier .mat exporte par export_tvlqr_for_bench.m")
    p.add_argument("--estimator", choices=["derivative", "ekf"],
                   default="ekf",
                   help="Kalman ETENDU (defaut) ou derivee filtree. Le filtre "
                        "lineaire d'estimator_v3 reste interdit : son rappel "
                        "change de signe au sommet du looping")
    p.add_argument("--trial", type=int, default=1)
    p.add_argument("--outdir", default=cfg.OUTDIR)
    p.add_argument("--dry-run", action="store_true",
                   help="charge et verifie le plan et la configuration, sans "
                        "toucher au materiel")
    args = p.parse_args()

    # Avant tout le reste, y compris --dry-run : un --dry-run qui ne verifie
    # pas la configuration donne une confiance immeritee.
    check_config_complete(cfg)

    plan = load_plan(args.plan)
    duration = plan["Tf"] + cfg.HOLD_AFTER_TF_S

    print("=" * 70)
    print("  T3 v3 -- looping, cerceau {}".format(cfg.HOOP))
    print("=" * 70)
    print("  plan       : {}".format(os.path.basename(args.plan)))
    print("  Tf         : {:.4f} s  (+ {:.1f} s de maintien)"
          .format(plan["Tf"], cfg.HOLD_AFTER_TF_S))
    print("  noeuds     : {} ({:.1f} Hz moyen, boucle a {:.0f} Hz)"
          .format(plan["t"].size, (plan["t"].size - 1) / plan["Tf"], cfg.LOOP_HZ))
    print("  psi_final  : {:+.1f} deg  (convention BANC, apres conversion "
          "dans l'export)".format(math.degrees(plan["x_ref"][-1, 2])))
    print("  max |u_ff| : {:.2f} rad/s^2   (U_SLIP_MAX = {:.0f})"
          .format(np.max(np.abs(plan["u_ff"])), cfg.U_SLIP_MAX))
    print("  max |th_d| : {:.2f} rad/s     (THETA_DOT_MAX = {:.0f})"
          .format(np.max(np.abs(plan["x_ref"][:, 1])), cfg.THETA_DOT_MAX))
    print("  max |psi_d|: {:.2f} rad/s     (file camera : {:.1f} mm a "
          "{:.0f} us de pose)"
          .format(np.max(np.abs(plan["x_ref"][:, 3])),
                  1e3 * cfg.R_EFF_M * np.max(np.abs(plan["x_ref"][:, 3]))
                  * cfg.EXPOSURE_US * 1e-6,
                  cfg.EXPOSURE_US))
    psid_max = float(np.max(np.abs(plan["x_ref"][:, 3])))
    print("  entre 2 img: {:.1f} mm parcourus ({:.2f} diametre de balle){}"
          .format(1e3 * cfg.R_EFF_M * psid_max / cfg.LOOP_HZ,
                  cfg.R_EFF_M * psid_max / cfg.LOOP_HZ / 0.025,
                  "" if getattr(cfg, "HOOP_RADIUS_PX", None) is None
                  else "  = {:.0f} px, fenetre {:.0f} px".format(
                      cfg.HOOP_RADIUS_PX * psid_max / cfg.LOOP_HZ,
                      cfg.TRACKING_WIN / 2.0)))
    print("  unwrap psi : {}".format(cfg.UNWRAP_PSI))
    print("  psi_dot    : {}".format(
        "Kalman etendu (sin psi) + compensation de retard"
        if args.estimator == "ekf"
        else "derivee filtree {:.0f} Hz".format(cfg.PSIDOT_FILTER_HZ)))
    print("  encodeur   : {:.0f} Hz (journal seul, hors loi de commande)"
          .format(cfg.LOOP_HZ / cfg.ENC_PERIOD))
    print("  theta_dot  : {}".format(
        "consigne integree (boucle interne supposee parfaite)"
        if cfg.TAU_VEL_S <= 0.0
        else "consigne retardee de tau = {:.1f} ms".format(1e3 * cfg.TAU_VEL_S)))

    problems = check_plan_against_config(plan)
    problems += geo.check_geometry(cfg)
    if problems:
        print("\n[!!] Le plan et la configuration du banc sont incompatibles :")
        for q in problems:
            print("     - {}".format(q))
        raise SystemExit("corriger t3_config_v3.py ou replanifier cote MATLAB")
    print("\n  plan compatible avec la configuration du banc")

    if args.dry_run:
        print("  --dry-run : rien n'a ete envoye au materiel")
        return

    bc.validate_config(cfg)

    stopper = bc.Stopper()
    picam2, det, roi, centre = bc.setup_camera(cfg)

    estimator, n_delay = build_estimator(args.estimator, cfg)
    if estimator is not None:
        print("[est] compensation de retard sur {} pas ({:.0f} ms)"
              .format(n_delay, 1e3 * n_delay * cfg.DT_NOMINAL))
        print("[est] sigma_psi = {:.3f} deg, q_accel = {:.1f}"
              .format(cfg.SIGMA_PSI_DEG, cfg.Q_ACCEL))
        print("[est] Q_ACCEL a ete regle au voisinage du fond ; pendant un "
              "looping l'erreur de")
        print("      modele dominante n'est plus la meme (force normale qui "
              "s'effondre au sommet).")

    odrv0 = ax = None
    try:
        odrv0, ax = bc.connect_and_prepare(cfg)

        # La balle demarre au repos au fond, comme dans le plan (x0 = 0).
        psi_offset, psi_init_meas, arm_stats = bc.arm(
            picam2, det, roi, centre, cfg, 0.0)

        log = bc.control_loop(
            ax, picam2, det, roi, centre, psi_offset,
            plan["K_traj"][0], trajectory_reference(plan), duration,
            stopper, cfg, estimator=estimator, tag="t3",
            gain_schedule=gain_schedule(plan))
        
        bc.save_and_summarise(
            log, cfg, args.outdir,
            "T3v3_{}_loop_{}_t{:02d}_{{stamp}}.npz".format(
                cfg.HOOP, args.estimator, args.trial),
            {"experiment": "T3", "version": 3, "hoop": cfg.HOOP,
             "trial": args.trial,
             "plan_file": os.path.basename(args.plan),
             "Tf": plan["Tf"],
             "psi_final_plan_deg": math.degrees(plan["x_ref"][-1, 2]),
             "arm_offset_deg": math.degrees(arm_stats["offset"]),
             "arm_noise_deg": math.degrees(arm_stats["sigma"]),
             "psi_init_measured_deg": psi_init_meas,
             "ref_type": "tvlqr_trajectory", "estimator": args.estimator,
             "estimator_n_delay": int(n_delay),
             "sigma_psi_deg": cfg.SIGMA_PSI_DEG, "q_accel": cfg.Q_ACCEL,
             "tracking_win": int(cfg.TRACKING_WIN),
             "hoop_radius_px": cfg.HOOP_RADIUS_PX,
             "roi_radius_px": getattr(cfg, "ROI_RADIUS_PX", None),
             "unwrap_psi": bool(cfg.UNWRAP_PSI),
             "psi_offset_rad": float(psi_offset),
             "roi": [int(v) for v in roi], },
            plan["K_traj"][0], "tvlqr")

    finally:
        bc.shutdown(ax, picam2, cfg)


if __name__ == "__main__":
    main()
