"""
e2_run.py -- Experience E2 : oscillations libres de la balle (sur le Pi).

Le hoop est bloque mecaniquement, l'ODrive reste en IDLE. Le script ne fait
que filmer la balle et enregistrer sa position angulaire psi image par image.
"""

import argparse
import json
import math
import os
import time
from datetime import datetime

import cv2
import numpy as np

from common.ball_detection import (BallDetectorAA4CC, open_camera,
                            DEFAULT_COLOR_COEFS)

# --- Réglages retenus à l'issue de la calibration caméra (C1-C5) ---------- #
SIZE          = (820, 616)
FPS           = 50
EXPOSURE_US   = 18000
GAIN          = 12.0
THRESHOLD     = 110        # valeur AA4CC d'origine, validée sur ce montage
DOWNSAMPLE    = 8
TRACKING_WIN  = 64
BALL_SIZE     = (0, 150)   # diamètre min/max [px]
ROI_MARGIN    = 0.05

# --- Paramètres Hough guidés par ta géométrie actuelle -------------------- #
EXPECTED_CENTER_PX    = (410.0, 310.0)
CENTER_TOLERANCE_PX   = 180.0
EXPECTED_OUTER_RADIUS_PX   = 280.0
OUTER_RADIUS_TOLERANCE_PX  = 45.0

HOUGH_DP        = 1.2
HOUGH_MIN_DIST  = 120
HOUGH_PARAM1    = 120
HOUGH_PARAM2    = 50
HOUGH_MIN_RADIUS = int(EXPECTED_OUTER_RADIUS_PX - OUTER_RADIUS_TOLERANCE_PX)
HOUGH_MAX_RADIUS = int(EXPECTED_OUTER_RADIUS_PX + OUTER_RADIUS_TOLERANCE_PX)


def detect_hough_circles(frame):
    """Détecte les cercles visibles dans une image complète."""
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    gray = cv2.GaussianBlur(gray, (9, 9), 2.0)

    circles = cv2.HoughCircles(
        gray,
        cv2.HOUGH_GRADIENT,
        dp=HOUGH_DP,
        minDist=HOUGH_MIN_DIST,
        param1=HOUGH_PARAM1,
        param2=HOUGH_PARAM2,
        minRadius=HOUGH_MIN_RADIUS,
        maxRadius=HOUGH_MAX_RADIUS,
    )

    if circles is None:
        return np.empty((0, 3), dtype=float)

    circles = np.asarray(circles[0], dtype=float)
    return circles[np.argsort(circles[:, 2])[::-1]]


def filter_by_center(circles, expected_center, tolerance):
    """Filtre les cercles dont le centre est proche de la position attendue."""
    if len(circles) == 0:
        return circles

    ex, ey = expected_center
    dist = np.hypot(circles[:, 0] - ex, circles[:, 1] - ey)
    return circles[dist <= tolerance]


def select_expected_circle(circles, expected_radius, tolerance):
    """Sélectionne le cercle dont le rayon est le plus proche de la valeur attendue."""
    if len(circles) == 0:
        return None

    dist = np.abs(circles[:, 2] - expected_radius)
    accepted = circles[dist <= tolerance]
    if len(accepted) == 0:
        return None

    idx = int(np.argmin(np.abs(accepted[:, 2] - expected_radius)))
    return accepted[idx]


def roi_from_circle(cx, cy, radius, margin=ROI_MARGIN, image_size=SIZE):
    """Construit une ROI rectangulaire autour du cercle détecté."""
    image_w, image_h = image_size
    half = radius * (1.0 + margin)

    x0 = max(0, int(round(cx - half)))
    y0 = max(0, int(round(cy - half)))
    x1 = min(image_w, int(round(cx + half)))
    y1 = min(image_h, int(round(cy + half)))

    return x0, y0, x1 - x0, y1 - y0


def run(duration_s, countdown_s, hoop, psi0_deg, trial, outdir):
    picam2 = open_camera(SIZE, FPS, EXPOSURE_US, GAIN)

    det = BallDetectorAA4CC(
        color_coefs=DEFAULT_COLOR_COEFS,
        threshold=THRESHOLD,
        downsample=DOWNSAMPLE,
        tracking_window=TRACKING_WIN,
        ball_size=BALL_SIZE,
    )

    # Échauffement
    for _ in range(10):
        picam2.capture_array()

    # --- Recadrage sur le hoop via Hough, une seule fois ---
    frame = picam2.capture_array()
    circles = detect_hough_circles(frame)

    if len(circles) == 0:
        picam2.stop()
        picam2.close()
        raise RuntimeError("E2: aucun cercle Hough détecté pour le hoop.")

    center_candidates = filter_by_center(
        circles,
        EXPECTED_CENTER_PX,
        CENTER_TOLERANCE_PX,
    )

    outer = select_expected_circle(
        center_candidates,
        EXPECTED_OUTER_RADIUS_PX,
        OUTER_RADIUS_TOLERANCE_PX,
    )

    if outer is None:
        picam2.stop()
        picam2.close()
        raise RuntimeError(
            "E2: aucun cercle compatible avec centre≈{} +/- {:.1f}px, "
            "rayon≈{:.1f} +/- {:.1f}px."
            .format(EXPECTED_CENTER_PX,
                    CENTER_TOLERANCE_PX,
                    EXPECTED_OUTER_RADIUS_PX,
                    OUTER_RADIUS_TOLERANCE_PX)
        )

    hcx, hcy, hradius = outer
    rx, ry, rw, rh = roi_from_circle(hcx, hcy, hradius)

    # Centre du hoop en coordonnées image complètes.
    cx, cy = rx + rw / 2.0, ry + rh / 2.0
    print("[e2] centre hoop (Hough) : ({:.0f}, {:.0f}), rayon~{:.0f}px"
          .format(cx, cy, hradius))
    print("[e2] ROI : x[{}:{}] y[{}:{}] ({} x {} px)"
          .format(rx, rx + rw, ry, ry + rh, rw, rh))

    # --- Compte à rebours ---
    print()
    for k in range(int(countdown_s), 0, -1):
        print("  lacher la balle dans {} ...".format(k))
        time.sleep(1.0)
    print("  >>> LACHER MAINTENANT <<<  ({:.0f} s d'enregistrement)\n"
          .format(duration_s))

    # --- Acquisition ---
    n_max = int(duration_s * FPS * 1.5) + 100     # marge confortable
    T   = np.empty(n_max)
    U   = np.full(n_max, np.nan)
    V   = np.full(n_max, np.nan)
    PSI = np.full(n_max, np.nan)
    OK  = np.zeros(n_max, dtype=bool)

    i = 0
    t0 = time.perf_counter()
    while True:
        t = time.perf_counter() - t0
        if t >= duration_s or i >= n_max:
            break

        frame = picam2.capture_array()
        loc = det.process_image(frame[ry:ry + rh, rx:rx + rw])

        T[i] = t
        if loc is not None:
            u, v = loc[0] + rx, loc[1] + ry
            U[i], V[i] = u, v
            PSI[i] = math.atan2(u - cx, v - cy)
            OK[i] = True
        i += 1

    picam2.stop()
    picam2.close()

    T, U, V, PSI, OK = T[:i], U[:i], V[:i], PSI[:i], OK[:i]

    # --- Rapport ---
    dt = np.diff(T)
    fps_eff = 1.0 / np.median(dt) if len(dt) else float('nan')
    rate = 100.0 * OK.mean() if OK.size else 0.0
    print("[e2] {} frames en {:.1f} s  ->  {:.1f} fps effectifs"
          .format(i, T[-1] if i else 0, fps_eff))
    print("[e2] detection : {:.1f} %".format(rate))
    if rate < 95:
        print("[!]  taux de detection faible : verifier l'eclairage / le seuil")

    if OK.sum() > 10:
        amp = np.degrees(np.nanmax(PSI[OK]) - np.nanmin(PSI[OK])) / 2
        print("[e2] amplitude crete observee : ~{:.1f} deg".format(amp))
        if amp < 10:
            print("[!]  amplitude faible : relacher la balle plus haut")

    # --- Sauvegarde ---
    os.makedirs(outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = "E2_{}_{:02.0f}deg_t{:02d}_{}.npz".format(hoop, psi0_deg, trial, stamp)
    path = os.path.join(outdir, name)

    meta = {
        "hoop": hoop, "psi0_nominal_deg": psi0_deg, "trial": trial,
        "timestamp": stamp, "fps_effective": float(fps_eff),
        "detection_rate_pct": float(rate),
        "hoop_center_px": [float(cx), float(cy)],
        "roi": [int(v) for v in (rx, ry, rw, rh)],
        "camera": {"size": list(SIZE), "fps": FPS,
                   "exposure_us": EXPOSURE_US, "gain": GAIN},
        "detector": {"threshold": THRESHOLD,
                     "color_coefs": list(DEFAULT_COLOR_COEFS),
                     "downsample": DOWNSAMPLE,
                     "tracking_window": TRACKING_WIN,
                     "ball_size": list(BALL_SIZE)},
        "hough": {
            "expected_center_px": list(EXPECTED_CENTER_PX),
            "expected_radius_px": float(EXPECTED_OUTER_RADIUS_PX),
            "selected_center_px": [float(hcx), float(hcy)],
            "selected_radius_px": float(hradius),
        },
    }

    np.savez_compressed(path, t=T, u=U, v=V, psi=PSI, ok=OK,
                        meta=json.dumps(meta))
    print("[out] {}".format(path))
    return path


def main():
    p = argparse.ArgumentParser(description="E2 : oscillations libres")
    p.add_argument("--duration", type=float, default=30.0,
                  help="duree d'enregistrement [s]")
    p.add_argument("--countdown", type=float, default=3.0)
    p.add_argument("--hoop", choices=["outer", "inner"], default="outer")
    p.add_argument("--psi0", type=float, default=50.0,
                  help="angle de lacher nominal [deg] (metadonnee)")
    p.add_argument("--trial", type=int, default=1)
    p.add_argument("--outdir", default="data/e2")
    args = p.parse_args()

    run(args.duration, args.countdown, args.hoop, args.psi0,
        args.trial, args.outdir)


if __name__ == "__main__":
    main()
    