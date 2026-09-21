#!/usr/bin/env python3
"""
sign_check_v3.py -- Determination experimentale du signe encoder vs caméra.

Principe
--------
Le modele lineaire autour de psi = 0 s'ecrit

    psi_ddot = A21 psi + A22 psi_dot + B2 u,      u = theta_ddot

Balle au repos au fond (psi = psi_dot = 0), on applique un echelon
d'acceleration u > 0 pendant une duree courte devant la periode propre
(1/f_n ~ 0.76 s). A l'instant initial, tous les termes sauf B2*u sont nuls,
donc

    psi_ddot(0+) = B2 * u

Le signe de la PREMIERE excursion de psi donne donc directement le signe de
B2. C'est une mesure de signe, pas d'amplitude : elle est insensible aux
erreurs sur A21, A22 et sur la valeur de |B2|.

Le test est fait en boucle OUVERTE : aucun gain n'intervient, donc aucune
possibilite qu'un signe de retour d'etat vienne masquer le resultat.

Securite
--------
u est limite a U_TEST et la duree a T_TEST : l'excursion de psi predite
reste sous ~10 deg et la vitesse du cerceau sous 2 rad/s. La rampe de
deceleration est appliquee dans tous les cas.

Usage
-----
    python3 sign_check_v3.py
"""

import math
import time

import numpy as np

from common import bench_common as bc
from tasks import t1_config as cfg

U_TEST = 4.0        # [rad/s^2] echelon d'acceleration
T_TEST = 0.30       # [s] duree de l'echelon
T_OBS = 0.60        # [s] duree totale d'observation


def main():
    print("=" * 70)
    print("  Verification du signe de B2 (boucle ouverte)")
    print("=" * 70)
    print("  B2 configure : {:+.3f}".format(cfg.B2))
    print("  echelon u = {:+.1f} rad/s^2 pendant {:.2f} s"
          .format(U_TEST, T_TEST))
    print("  excursion initiale de psi attendue : POSITIVE si B2 > 0,")
    print("                                       NEGATIVE si B2 < 0")

    stopper = bc.Stopper()
    picam2, det, roi, centre = bc.setup_camera(cfg)
    rx, ry, rw, rh = roi
    cx, cy = centre

    odrv0 = ax = None
    theta_dot_cmd = 0.0
    try:
        odrv0, ax = bc.connect_and_prepare(cfg)
        psi_offset, _, _ = bc.arm(picam2, det, roi, centre, cfg, 0.0)

        ts, psis = [], []
        t_start = time.perf_counter()
        t_prev = t_start
        while True:
            frame = picam2.capture_array()
            t = time.perf_counter()
            dt = t - t_prev
            t_prev = t
            t_rel = t - t_start
            if t_rel > T_OBS or stopper.stop:
                break

            u = U_TEST if t_rel < T_TEST else 0.0
            theta_dot_cmd += u * dt
            ax.controller.input_vel = theta_dot_cmd / (2 * math.pi)

            loc = det.process_image(frame[ry:ry + rh, rx:rx + rw])
            if loc is None:
                continue
            u_px, v_px = loc[0] + rx, loc[1] + ry
            psis.append(math.atan2(u_px - cx, v_px - cy) - psi_offset)
            ts.append(t_rel)

        ts, psis = np.asarray(ts), np.asarray(psis)
        if ts.size < 10:
            raise RuntimeError("pas assez de detections.")

        # On regarde l'excursion pendant l'echelon seulement : au-dela, le
        # rappel gravitaire ramene la balle et le signe s'inverse.
        m = ts < T_TEST
        peak = psis[m][np.argmax(np.abs(psis[m]))]
        print("\n  excursion max pendant l'echelon : {:+.2f} deg"
              .format(math.degrees(peak)))
        b2_sign = 1.0 if peak > 0 else -1.0
        print("  => signe de B2 mesure : {:+.0f}".format(b2_sign))
        if b2_sign * cfg.B2 > 0:
            print("  => COHERENT avec la configuration ({:+.3f}). Rien a "
                  "changer.".format(cfg.B2))
        else:
            print("  => INCOHERENT avec la configuration ({:+.3f})."
                  .format(cfg.B2))
            print("     Mettre B2 = {:+.3f} dans t1_config_v3.py AVANT de "
                  "lancer T2.".format(-cfg.B2))
        if abs(math.degrees(peak)) < 0.8:
            print("\n  [!] excursion faible devant le bruit de mesure "
                  "(sigma_psi = {:.2f} deg) : refaire avec U_TEST plus grand."
                  .format(cfg.SIGMA_PSI_DEG))

    finally:
        if ax is not None:
            bc.ramp_down(ax, theta_dot_cmd, cfg)
        bc.shutdown(ax, picam2)


if __name__ == "__main__":
    main()
