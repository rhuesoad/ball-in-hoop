#!/usr/bin/env python3

import argparse
import math
import numpy as np
from common import bench_common as bc
from tasks import t1_config as cfg


def build_estimator(name, cfg):
    if name == "derivative":
        return None, 0
    if name == "kalman":
        from common import estimator as est_mod
    elif name == "ekf":
        from common import ekf as est_mod
    else:
        raise SystemExit("estimateur inconnu : {}".format(name))
    return est_mod.build_estimator(cfg, cfg.DT_NOMINAL)


def constant_reference(psi_ref):
    """Reference T1 : position de consigne fixe, aucun feedforward."""
    x_ref = np.array([0.0, 0.0, psi_ref, 0.0])

    def ref(t):
        return x_ref, cfg.U_EQ
    return ref


def main():
    # ARGUMENTS 
    p = argparse.ArgumentParser(
        description="T1 : stabilisation de la balle")
    p.add_argument("--psi-init", type=float, default=0.0,
                   help="position de lacher de la balle [deg]")
    p.add_argument("--psi-ref", type=float, default=0.0,
                   help="consigne de position [deg]")
    p.add_argument("--gain", default=cfg.DEFAULT_GAIN,
                   choices=sorted(cfg.GAINS),
                   help="jeu de gains LQR")
    p.add_argument("--estimator", choices=["derivative", "kalman", "ekf"],
                   default="derivative",
                   help="estimateur de psi_dot (defaut : derivee filtree)")
    p.add_argument("--duration", type=float, default=30.0)
    p.add_argument("--trial", type=int, default=1)
    p.add_argument("--outdir", default=cfg.OUTDIR)
    args = p.parse_args()

    if args.duration > cfg.MAX_DURATION_S:
        raise SystemExit("duree superieure au plafond de {:.0f} s"
                         .format(cfg.MAX_DURATION_S))

    # Validate the configuration 
    bc.validate_config(cfg)

    # Pick the right set of gains 
    K = bc.pick_gain(cfg, args.gain)

    psi_ref = math.radians(args.psi_ref)

    print("=" * 70)
    print("  T1 v3 -- stabilisation, cerceau {}".format(cfg.HOOP))
    print("=" * 70)
    print("  gains ({})  : K = [{:.4f}, {:.4f}, {:.4f}, {:.4f}]"
          .format(args.gain, *K))
    print("  psi_init = {:+.1f} deg,  psi_ref = {:+.1f} deg"
          .format(args.psi_init, args.psi_ref))
    print("  theta    : {}".format(
        "encodeur" if cfg.USE_MEASURED_THETA else "integration de u"))
    print("  psi_dot  : {}".format(
        {"kalman": "Kalman lineaire + compensation de retard",
         "ekf": "Kalman etendu (sin psi) + compensation de retard"}
        .get(args.estimator,
             "derivee filtree {:.0f} Hz".format(cfg.PSIDOT_FILTER_HZ))))
    if cfg.GAIN_START < 1.0:
        print("  gain     : {:.2f} -> 1.00 en {:.2f} s"
              .format(cfg.GAIN_START, cfg.GAIN_RAMP_S))
    else:
        print("  gain     : nominal des le premier pas")
    print("  ODrive   : consigne {:.0f} Hz, encodeur {:.0f} Hz, "
          "erreurs {:.0f} Hz".format(cfg.LOOP_HZ / cfg.COM_PERIOD,
                                     cfg.LOOP_HZ / cfg.ENC_PERIOD,
                                     cfg.LOOP_HZ / cfg.ERR_PERIOD))
    print("  modele   : f_n = {:.4f} Hz, zeta = {:.3f}, B2 = {:+.3f}"
          .format(cfg.F_N_HZ, cfg.ZETA, cfg.B2))

    # Verification de stabilite AVANT tout mouvement. Elle ne coute rien et
    # signale immediatement une incoherence de signe entre K et B2.
    A = np.array([[0, 1, 0, 0], [0, 0, 0, 0],
                  [0, 0, 0, 1], [0, 0, cfg.A21, cfg.A22]], float)
    B = np.array([0.0, 1.0, 0.0, cfg.B2])
    lam = np.linalg.eigvals(A - np.outer(B, K))
    print("  poles boucle fermee (modele) : max Re = {:+.3f}"
          .format(float(lam.real.max())))
    if lam.real.max() > 0:
        print("\n[!!] Le modele predit une boucle fermee INSTABLE avec ce jeu")
        print("     de gains et ce signe de B2. Ne pas lancer sans avoir")
        print("     verifie le signe de B2 sur le banc (sign_check_v3.py).")
        if input("     Continuer quand meme ? [oui/non] ").strip() != "oui":
            raise SystemExit("annule")

    # Constante de temps du mode le plus lent (c'est theta qui traine, pas
    # psi). La "seconde moitie" utilisee par save_and_summarise n'a de sens
    # que si le transitoire y est termine.
    tau_slow = -1.0 / float(lam.real.max())
    print("  mode le plus lent : tau = {:.1f} s".format(tau_slow))
    if args.duration < 5.0 * tau_slow:
        print("\n  [!] duree de {:.0f} s pour tau = {:.1f} s : le transitoire de"
              .format(args.duration, tau_slow))
        print("      theta n'est pas termine dans la seconde moitie de l'essai,")
        print("      que les statistiques de fin utilisent. Duree conseillee :")
        print("      {:.0f} s au moins.".format(math.ceil(5.0 * tau_slow)))

    # Effort demande par la condition initiale, AVANT tout mouvement. A
    # l'instant du lacher psi = psi_init et tout le reste est nul, donc
    # u(0+) = -K3 * psi_init exactement. Si ce seul terme sature deja, le
    # retour d'etat ne fait pas ce pour quoi il a ete synthetise et l'essai
    # ne compare plus rien a la simulation.
    u0 = abs(K[2]) * abs(math.radians(args.psi_init))
    print("  |u| au lacher : {:.1f} rad/s^2  (borne de glissement {:.0f})"
          .format(u0, cfg.U_SLIP_MAX))

    
    stopper = bc.Stopper()
    picam2, det, roi, centre = bc.setup_camera(cfg)         # Camera calibration, verification that the centre is aligned & hoop has the right scale.

    estimator, n_delay = build_estimator(args.estimator, cfg)
    if estimator is not None:
        print("[est] compensation de retard sur {} pas ({:.0f} ms)"
              .format(n_delay, 1e3 * n_delay * cfg.DT_NOMINAL))
        print("[est] sigma_psi = {:.3f} deg, q_accel = {:.1f}"
              .format(cfg.SIGMA_PSI_DEG, cfg.Q_ACCEL))

    odrv0 = ax = None
    try:
        odrv0, ax = bc.connect_and_prepare(cfg)
        psi_offset, psi_init_meas, arm_stats = bc.arm(
            picam2, det, roi, centre, cfg, args.psi_init)

        log = bc.control_loop(
            ax, picam2, det, roi, centre, psi_offset,
            K, constant_reference(psi_ref), args.duration, stopper, cfg,
            estimator=estimator, tag="t1")

        bc.save_and_summarise(
            log, cfg, args.outdir,
            "T1v3_{}_init{:+03.0f}_ref{:+03.0f}_{}_{}_t{:02d}_{{stamp}}.npz"
            .format(cfg.HOOP, args.psi_init, args.psi_ref, args.gain,
                    args.estimator, args.trial),
            {"experiment": "T1", "version": 3, "hoop": cfg.HOOP,
             "trial": args.trial,
             "psi_init_deg": args.psi_init,
             "psi_init_measured_deg": psi_init_meas,
             "arm_offset_deg": math.degrees(arm_stats["offset"]),
             "arm_osc_deg": math.degrees(arm_stats["osc_amp"]),
             "arm_noise_deg": math.degrees(arm_stats["sigma"]),
             "psi_ref_deg": args.psi_ref,
             "ref_type": "constant", "estimator": args.estimator,
             "estimator_n_delay": int(n_delay),
             "sigma_psi_deg": cfg.SIGMA_PSI_DEG, "q_accel": cfg.Q_ACCEL,
             "psi_offset_rad": float(psi_offset),
             "roi": [int(v) for v in roi]},
            K, args.gain)

    finally:
        bc.shutdown(ax, picam2, cfg)


if __name__ == "__main__":
    main()
