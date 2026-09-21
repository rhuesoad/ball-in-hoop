#!/usr/bin/env python3
"""
inner_loop_step.py -- Caracterisation de la boucle de vitesse de l'ODrive.

L'architecture en cascade suppose que la boucle interne de vitesse est
beaucoup plus rapide que la boucle externe. Ce script mesure a quel point
c'est vrai : il envoie des echelons de vitesse et enregistre la reponse
aussi vite que l'USB le permet.

Ce qu'on cherche
----------------
  - le temps de montee a 90 % : c'est lui qui dit si la separation
    d'echelles de temps tient. La balle oscille avec une periode de 760 ms ;
    il faudrait idealement moins de 40 ms ici.
  - la forme de la reponse : un premier ordre (montee douce, sans
    depassement) se compense en accelerant les gains ; une reponse avec
    depassement ou plateau signale autre chose (limite de courant,
    integrateur sature).
  - les gains actuels vel_gain et vel_integrator_gain : s'ils sont aux
    valeurs par defaut, ils n'ont jamais ete accordes pour cette inertie.

Le cerceau TOURNE pendant ce test. Degager la zone. La balle peut rester
en place, elle ne sera pas propulsee bien loin a 1 tr/s.

Usage
-----
    python3 inner_loop_step.py
    python3 inner_loop_step.py --vel 2.0 --repeats 5
"""

import argparse
import json
import math
import os
import time
from datetime import datetime

import numpy as np
import odrive
from odrive.enums import *


def connect():
    print("[odrv] connexion ...")
    d = odrive.find_any()
    ax = d.axis0
    print("[odrv] connecte (serial {}, fw {}.{}.{})".format(
        d.serial_number, d.fw_version_major,
        d.fw_version_minor, d.fw_version_revision))

    if not d.config.brake_resistor0.enable:
        raise RuntimeError("resistance de freinage desactivee")

    if ax.current_state != AXIS_STATE_CLOSED_LOOP_CONTROL:
        print("[odrv] calibration de l'offset encodeur ...")
        ax.requested_state = AXIS_STATE_ENCODER_OFFSET_CALIBRATION
        while ax.current_state != AXIS_STATE_IDLE:
            time.sleep(0.1)
        if ax.active_errors or ax.disarm_reason:
            raise RuntimeError("calibration echouee : active=0x{:X}"
                               .format(ax.active_errors))
        ax.requested_state = AXIS_STATE_CLOSED_LOOP_CONTROL
        time.sleep(0.5)

    ax.controller.config.control_mode = CONTROL_MODE_VELOCITY_CONTROL
    ax.controller.config.input_mode = INPUT_MODE_PASSTHROUGH
    ax.controller.input_vel = 0.0
    return d, ax


def one_step(ax, vel_tps, hold_s, settle_s):
    """Un echelon, echantillonne aussi vite que l'USB le permet."""
    ax.controller.input_vel = 0.0
    time.sleep(settle_s)

    T, V = [], []
    t0 = time.perf_counter()
    ax.controller.input_vel = vel_tps
    while True:
        t = time.perf_counter() - t0
        if t >= hold_s:
            break
        T.append(t)
        V.append(ax.vel_estimate)

    ax.controller.input_vel = 0.0
    time.sleep(settle_s)
    return np.array(T), np.array(V)


def analyse(T, V, vel_tps):
    """Temps de montee, depassement, constante de temps equivalente."""
    out = {"fs_hz": 1.0 / np.median(np.diff(T)), "n": T.size}

    def crossing(frac):
        idx = np.argmax(V >= frac * vel_tps)
        if V[idx] < frac * vel_tps:
            return float("nan")
        if idx == 0:
            return T[0]
        # interpolation lineaire entre les deux echantillons encadrants
        v0, v1 = V[idx - 1], V[idx]
        return T[idx - 1] + (T[idx] - T[idx - 1]) * \
            (frac * vel_tps - v0) / (v1 - v0)

    out["t10"] = crossing(0.10)
    out["t63"] = crossing(0.632)      # constante de temps d'un premier ordre
    out["t90"] = crossing(0.90)
    out["overshoot_pct"] = 100.0 * (V.max() - vel_tps) / vel_tps
    tail = V[T > 0.8 * T[-1]]
    out["steady_gain"] = float(tail.mean() / vel_tps) if tail.size else float("nan")
    return out


def main():
    p = argparse.ArgumentParser(
        description="Reponse indicielle de la boucle de vitesse ODrive")
    p.add_argument("--vel", type=float, default=1.0,
                   help="amplitude de l'echelon [tr/s] (defaut 1.0)")
    p.add_argument("--hold", type=float, default=0.6,
                   help="duree d'enregistrement apres l'echelon [s]")
    p.add_argument("--settle", type=float, default=1.0,
                   help="pause a vitesse nulle entre les essais [s]")
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--outdir", default="data/inner_loop")
    args = p.parse_args()

    d, ax = connect()

    print("\n" + "=" * 62)
    print("  Reglages actuels de la boucle de vitesse")
    print("=" * 62)
    cc = ax.controller.config
    for name in ("vel_gain", "vel_integrator_gain", "vel_limit",
                 "pos_gain", "vel_integrator_limit"):
        print("  {:<24s} = {}".format(name, getattr(cc, name, "ABSENT")))
    print("  {:<24s} = {}".format("current_soft_max",
                                  ax.config.motor.current_soft_max))
    print("  {:<24s} = {:.2f} V".format("vbus", d.vbus_voltage))

    input("\nZone degagee, le cerceau va tourner. ENTREE pour lancer ...")

    runs, stats = [], []
    try:
        for i in range(args.repeats):
            T, V = one_step(ax, args.vel, args.hold, args.settle)
            s = analyse(T, V, args.vel)
            runs.append((T, V))
            stats.append(s)
            print("  essai {}: fs={:5.0f} Hz  t10={:5.1f} ms  t63={:5.1f} ms  "
                  "t90={:5.1f} ms  depassement={:+5.1f} %  gain={:.3f}"
                  .format(i + 1, s["fs_hz"], 1e3 * s["t10"], 1e3 * s["t63"],
                          1e3 * s["t90"], s["overshoot_pct"],
                          s["steady_gain"]))
            if ax.active_errors or ax.disarm_reason:
                print("[!] erreur ODrive : active=0x{:X}, disarm=0x{:X}"
                      .format(ax.active_errors, ax.disarm_reason))
                break
    finally:
        ax.controller.input_vel = 0.0
        time.sleep(0.3)
        ax.requested_state = AXIS_STATE_IDLE
        print("\n[odrv] IDLE")

    if not stats:
        return

    t90 = np.array([s["t90"] for s in stats])
    t63 = np.array([s["t63"] for s in stats])

    print("\n" + "=" * 62)
    print("  SYNTHESE")
    print("=" * 62)
    print("  t90 = {:.1f} +/- {:.1f} ms".format(1e3 * t90.mean(),
                                                1e3 * t90.std(ddof=1)
                                                if t90.size > 1 else 0.0))
    print("  t63 = {:.1f} ms  -> constante de temps equivalente"
          .format(1e3 * np.nanmean(t63)))
    print("  cadence d'echantillonnage : {:.0f} Hz"
          .format(np.mean([s["fs_hz"] for s in stats])))

    # Le pendule oscille a 8.31 rad/s (periode 756 ms). La separation
    # d'echelles est confortable si la boucle interne est au moins 20 fois
    # plus rapide, soit t90 <= 38 ms environ.
    tau = np.nanmean(t63)
    ratio = 0.756 / max(1e-6, np.nanmean(t90))
    print("\n  periode propre de la balle : 756 ms")
    print("  rapport periode / t90      : {:.1f}".format(ratio))
    if ratio >= 20:
        print("  -> separation d'echelles confortable, la boucle interne")
        print("     peut etre negligee dans la conception du LQR.")
    elif ratio >= 8:
        print("  -> separation marginale : la boucle interne consomme de la")
        print("     marge de phase. Accorder vel_gain, ou l'integrer au modele.")
    else:
        print("  -> separation INSUFFISANTE : la boucle interne fait partie")
        print("     de la dynamique a controler. Aucun jeu de gains concu en")
        print("     la negligeant ne tiendra.")
    print("\n  dephasage apporte a 8.31 rad/s (1er ordre, tau={:.0f} ms) : "
          "{:.0f} deg".format(1e3 * tau, math.degrees(math.atan(8.31 * tau))))

    os.makedirs(args.outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(args.outdir, "inner_loop_step_{}.npz".format(stamp))
    np.savez_compressed(
        path,
        meta=json.dumps({
            "vel_step_tps": args.vel, "repeats": len(runs),
            "vel_gain": cc.vel_gain,
            "vel_integrator_gain": cc.vel_integrator_gain,
            "stats": stats,
        }),
        **{"t{}".format(i): r[0] for i, r in enumerate(runs)},
        **{"v{}".format(i): r[1] for i, r in enumerate(runs)})
    print("\n[out] {}".format(path))


if __name__ == "__main__":
    main()
