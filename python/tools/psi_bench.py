#!/usr/bin/env python3
"""
psi_bench_v3.py -- Caracterisation de la mesure de psi, balle immobile.

Rien ne bouge, aucun risque. A lancer chaque fois qu'un doute existe 
sur la vision, et une fois avant chaque campagne pour avoir la valeur
de sigma_psi du jour a citer dans le memoire.

Ce qui est mesuré
-----------------

La balle est immobile au fond du cerceau : la physique impose psi = 0 et
psi_dot = 0. Tout ce qui est observe est donc soit un biais, du bruit ou 
une osocillation résiduelle. Le script sépare les 3 via une régression 
linéaire de psi(t) sur [sin(w_n t), cos(w_n t), 1], puis analyse le résidu.

Sorties
-------
    offset          biais systematique  -> qualité de HOOP_CENTRE_PX
    oscillation     amplitude a f_n     -> la balle bouge-t-elle vraiment
    sigma (MAD)     bruit robuste       -> a comparer a SIGMA_PSI_DEG (E3)
    sigma (std)     bruit non robuste   -> l'ecart avec la MAD chiffre le
                                           poids des detections aberrantes
    p2p, aberrants  amplitude et taux des valeurs extremes
    derive lineaire pente sur la duree  -> banc ou camera qui bougeZZZ
"""

import argparse
import math
import os
import time
from datetime import datetime

import numpy as np

from common import bench_common as bc
from tasks import t1_config as cfg


def main():
    p = argparse.ArgumentParser(description="Caracterisation de la mesure de psi (balle immobile)")
    p.add_argument("--duration", type=float, default=10.0)
    p.add_argument("--outdir", default="data/diag")
    args = p.parse_args()

    print("=" * 70)
    print("  psi_bench -- qualite de la mesure de psi (aucun mouvement)")
    print("=" * 70)
    print("  centre utilisé : ({:.1f}, {:.1f}) px  [{}]".format(
        cfg.HOOP_CENTRE_PX[0], cfg.HOOP_CENTRE_PX[1],
        "Hough" if cfg.USE_HOUGH_CENTRE else "Taubin/E3"))
    print("  reference E3   : sigma_psi = {:.2f} deg".format(cfg.SIGMA_PSI_DEG))

    picam2, det, roi, centre = bc.setup_camera(cfg)
    rx, ry, rw, rh = roi
    cx, cy = centre

    print("\n  Posez la balle AU FOND du cerceau et laissez-la s'immobiliser")
    print("  completement (compter ~15 s apres le dernier contact).")
    input("  ENTREE pour lancer l'acquisition de {:.0f} s ...".format(
        args.duration))

    ts, psis, us, vs = [], [], [], []
    n_miss = 0
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < args.duration:
        frame = picam2.capture_array()
        t = time.perf_counter() - t0
        loc = det.process_image(frame[ry:ry + rh, rx:rx + rw])
        if loc is None:
            n_miss += 1
            continue
        u_px, v_px = loc[0] + rx, loc[1] + ry
        ts.append(t)
        psis.append(math.atan2(u_px - cx, v_px - cy))
        us.append(u_px)
        vs.append(v_px)

    try:
        picam2.stop()
        picam2.close()
    except Exception:
        pass

    t = np.asarray(ts)
    psi = np.asarray(psis)
    u = np.asarray(us)
    v = np.asarray(vs)
    if t.size < 50:
        raise SystemExit("trop peu de detections ({}).".format(t.size))

    st = bc.psi_statistics(t, psi, cfg.F_N_HZ)

    # Residu apres retrait de l'oscillation et du biais
    w = 2 * math.pi * cfg.F_N_HZ
    M = np.column_stack([np.sin(w * t), np.cos(w * t), np.ones(t.size)])
    coef, *_ = np.linalg.lstsq(M, psi, rcond=None)
    resid = psi - M @ coef

    # Derive lineaire, testee sur le residu
    Md = np.column_stack([t, np.ones(t.size)])
    cd, *_ = np.linalg.lstsq(Md, resid, rcond=None)
    slope_deg_per_min = math.degrees(cd[0]) * 60.0

    d = np.degrees
    print("\n" + "=" * 70)
    print("  RESULTATS  ({} images, {} non detectees, {:.1f} % de perte)"
          .format(t.size, n_miss, 100.0 * n_miss / (t.size + n_miss)))
    print("=" * 70)
    print("  cadence effective   : {:6.1f} Hz".format(t.size / args.duration))
    print("\n  -- biais --")
    print("  offset (regression) : {:+7.3f} deg   (limite armement {:.1f})"
          .format(d(st["offset"]), cfg.ARM_PSI_MAX_DEG))
    print("  mediane brute       : {:+7.3f} deg".format(d(st["median"])))
    print("  centre implique     : cx = {:.1f} px  (configure {:.1f})"
          .format(np.median(u) - np.median(v - cy) * math.tan(st["offset"]),
                  cx))
    print("\n  -- mouvement residuel --")
    print("  oscillation a {:.2f} Hz: {:7.3f} deg   (limite armement {:.1f})"
          .format(cfg.F_N_HZ, d(st["osc_amp"]), cfg.ARM_OSC_MAX_DEG))
    print("\n  -- bruit de detection --")
    print("  sigma robuste (MAD) : {:7.3f} deg   (limite armement {:.1f})"
          .format(d(st["sigma"]), cfg.ARM_NOISE_MAX_DEG))
    print("  sigma brut (std)    : {:7.3f} deg".format(d(np.std(resid))))
    print("  crete-a-crete       : {:7.3f} deg".format(d(st["p2p"])))
    print("  aberrants (>5 sigma): {:6.1f} %".format(100 * st["outlier_frac"]))
    print("  en pixels           : sigma_u = {:.2f} px, sigma_v = {:.2f} px"
          .format(float(np.std(u)), float(np.std(v))))
    print("\n  -- derive --")
    print("  pente               : {:+7.3f} deg/min".format(slope_deg_per_min))

    print("\n" + "-" * 70)
    ratio = np.std(resid) / st["sigma"] if st["sigma"] > 0 else 1.0
    if d(st["sigma"]) > cfg.ARM_NOISE_MAX_DEG:
        print("  [!!] bruit au-dessus du seuil d'armement : la mesure de psi")
        print("       est inexploitable en l'etat. Regler la vision AVANT")
        print("       toute campagne (eclairage, reflets, THRESHOLD,")
        print("       BALL_SIZE, DOWNSAMPLE).")
    elif d(st["sigma"]) > 3 * cfg.SIGMA_PSI_DEG:
        print("  [!]  bruit {:.1f}x superieur a la valeur E3 ({:.2f} deg) :"
              .format(d(st["sigma"]) / cfg.SIGMA_PSI_DEG, cfg.SIGMA_PSI_DEG))
        print("       exploitable, mais SIGMA_PSI_DEG est a remettre a jour")
        print("       dans la configuration et dans le memoire.")
    else:
        print("  vision conforme a la campagne E3.")
    if ratio > 2.0:
        print("  [!]  sigma brut {:.1f}x sigma robuste : detections aberrantes"
              .format(ratio))
        print("       isolees. Chercher des reflets sur le cerceau ou une")
        print("       seconde zone rouge dans le champ.")
    if abs(slope_deg_per_min) > 0.5:
        print("  [!]  derive de {:+.2f} deg/min : montage mecanique ou camera"
              .format(slope_deg_per_min))
        print("       qui bouge, ou eclairage qui change.")

    os.makedirs(args.outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(args.outdir, "psi_bench_{}.npz".format(stamp))
    np.savez_compressed(path, t=t, psi=psi, u=u, v=v,
                        centre=np.array(centre), n_miss=n_miss)
    print("\n[out] {}".format(path))


if __name__ == "__main__":
    main()
