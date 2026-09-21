#!/usr/bin/env python3
"""
identify_motor.py — Acquisition pour l'identification J, b, tau_c.

Experience 1 a et b: moteur seul et moteur + cerceaux 

==> Applique une série d'échelons de couple depuis 0 dans les 2 sens et enregistre pos/vel/Iq bruts. 

==> Identification du modèle sans la bille --> Trouver J, b, tau_c par régression linéaire sur les échelons. 
    J*theta_ddot + b*theta_dot + tau_c*sign(theta_dot) = tau

==> Analyse faite dans analyze_identification.py

Usage :
    python3 identify_motor.py --probe            # balayage grossier
    python3 identify_motor.py --run E1a          # séquence complète, moteur nu
    python3 identify_motor.py --run E1b          # séquence complète, avec cerceaux

Convention d'unités :
    - Vitesses en tours et tours/s (ODrive), conversion en rad/s à l'enregistrement
    - Couples en N.m (ODrive) 
"""

import argparse         # Argparse pour arguments de ligne de commande
import time
import os
import sys
import numpy as np
import odrive
from odrive.enums import (
    AXIS_STATE_CLOSED_LOOP_CONTROL,
    AXIS_STATE_IDLE,
    CONTROL_MODE_TORQUE_CONTROL,
    INPUT_MODE_PASSTHROUGH,
)

TWO_PI = 2.0 * np.pi

# ----------------------------------------------------------------------
# Paramètres des expériences : à régler selon le moteur et la charge 
# ----------------------------------------------------------------------
# Amplitudes de couple [N.m], sens positif. Tests effectués sur - tau aussi. 
TAU_LIST = [0.02, 0.04, 0.06, 0.08, 0.10]

SETTLE_HORIZON = 6.0      # durée de log par échelon [s]
PAUSE          = 2.0      # roue libre entre deux échelons [s]
VEL_LIMIT      = 40.0     # vitesse limite de l'ODrive [tr/s] (pour éviter de dépasser la limite de l'ODrive)
LOG_DT_TARGET  = 0.002    # période de polling visée [s] (500 Hz)
ZERO_VEL_TOL   = 0.05     # |theta_dot| considéré "à l'arrêt" [rad/s]

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs_id")


# ----------------------------------------------------------------------
# Connexion / configuration
# ----------------------------------------------------------------------
def connect():
    print("Connexion à l'ODrive...", flush=True)
    odrv = odrive.find_any()
    ax = odrv.axis0
    if ax.encoder.is_ready is False:
        print("  ATTENTION : l'encodeur n'est pas prêt (calibration faite ?)")
    print(f"  Vbus = {odrv.vbus_voltage:.1f} V   "
          f"Kt = {ax.motor.config.torque_constant:.5f} N.m/A")
    return odrv, ax


def enter_torque_mode(ax):
    ax.controller.config.control_mode = CONTROL_MODE_TORQUE_CONTROL
    ax.controller.config.input_mode = INPUT_MODE_PASSTHROUGH
    ax.controller.config.vel_limit = VEL_LIMIT
    ax.controller.input_torque = 0.0
    ax.requested_state = AXIS_STATE_CLOSED_LOOP_CONTROL
    time.sleep(0.2)         # temps pour que ça se connecte bien


def stop(ax):
    ax.controller.input_torque = 0.0
    ax.requested_state = AXIS_STATE_IDLE


# ----------------------------------------------------------------------
# Mesures
# ----------------------------------------------------------------------
def read_sample(ax):
    """Echantillon brut, pos/vel convertis en rad, rad/s."""
    pos = ax.encoder.pos_estimate * TWO_PI      # rad
    vel = ax.encoder.vel_estimate * TWO_PI      # rad/s (FILTRÉ par l'observateur)
    iq  = ax.motor.current_control.Iq_measured  # A
    return pos, vel, iq


def wait_until_stopped(ax, timeout=5.0):
    """Roue libre jusqu'à l'arrêt effectif."""
    ax.controller.input_torque = 0.0
    t0 = time.time()
    while time.time() - t0 < timeout:
        _, vel, _ = read_sample(ax)
        if abs(vel) < ZERO_VEL_TOL:
            return True
        time.sleep(0.02)
    return False


def log_step(ax, tau, horizon):
    """
    Applique un échelon de couple et logge jusqu'à horizon.
    Retour : tableaux numpy (t, pos, vel, iq).
    """
    t_buf, p_buf, v_buf, i_buf = [], [], [], []
    ax.controller.input_torque = float(tau)
    t_start = time.perf_counter()
    while True:
        t = time.perf_counter() - t_start
        if t > horizon:
            break
        pos, vel, iq = read_sample(ax)
        t_buf.append(t); p_buf.append(pos); v_buf.append(vel); i_buf.append(iq)
        # garde-fou : la vel_limit de l'ODrive coupe déjà, mais on double
        if abs(vel) > 0.95 * VEL_LIMIT * TWO_PI:
            print("    !! vel_limit approchée, arrêt de l'échelon")
            break
        # cadence : busy-loop léger sur la cible
        next_t = t + LOG_DT_TARGET
        while time.perf_counter() - t_start < next_t:
            pass
    ax.controller.input_torque = 0.0
    return {
        "t":   np.asarray(t_buf),
        "pos": np.asarray(p_buf),
        "vel": np.asarray(v_buf),
        "iq":  np.asarray(i_buf),
        "tau_cmd": float(tau),
    }


def actual_rate(step):
    """Cadence effective atteinte [Hz], pour vérifier la faisabilité de T."""
    dt = np.diff(step["t"])
    return 1.0 / np.median(dt) if len(dt) else np.nan


# ----------------------------------------------------------------------
# Mode PROBE : ordre de grandeur de tau_c et de la vitesse atteinte
# ----------------------------------------------------------------------
def probe(ax):
    print("\n=== DECOLLAGE: balayage jusqu'au démarrage ===")
    print("But : trouver ~tau_c (décollage) et calibrer TAU_LIST.\n")

    # 1) Décollage : rampe lente jusqu'au premier mouvement du moteur.
    print("[1] Couple de décollage (breakaway)...")
    tau = 0.0
    dtau = 0.001
    wait_until_stopped(ax)              # moteur bien à l'arrêt 
    pos0, _, _ = read_sample(ax)        # position de référence

    ax.controller.input_torque = 0.0 
    t0 = time.perf_counter()
    tau_break = None

    # On monte le couple petit à petit jusqu'au premier mouvement.
    while tau < 0.3:
        tau += dtau
        ax.controller.input_torque = tau

        time.sleep(0.05)

        pos, vel, _ = read_sample(ax)

        # Si mouvement, on arrête. 
        if abs(pos - pos0) > 0.02 or abs(vel) > 0.1:   # premier mouvement
            tau_break = tau
            break
    
    ax.controller.input_torque = 0.0
    if tau_break:
        print(f"    tau_s (décollage) ~= {tau_break:.4f} N.m")
        print(f"    -> borne basse TAU_LIST conseillée ~ {2*tau_break:.3f} N.m")
    else:
        print("    pas de décollage détecté sous 0.3 N.m :(")
    wait_until_stopped(ax)


    # 2) Quels couples doit-on parcourir pour déterminer les paramètres dans la plage qui nous intéresse? 
    print("\n[2] Vitesses de régime pour couples test ")
    for tau_test in (0.05, 0.10):
        wait_until_stopped(ax)
        s = log_step(ax, tau_test, horizon=4.0)
        v_inf = np.median(s["vel"][-max(5, len(s["vel"])//10):])
        print(f"    tau={tau_test:.3f} N.m -> theta_dot_inf ~= {v_inf:6.2f} rad/s"
              f"   (cadence {actual_rate(s):.0f} Hz)")
    wait_until_stopped(ax)
    print("\n Fini! \n")


# ----------------------------------------------------------------------
# Mode RUN : séquence complète
# ----------------------------------------------------------------------
def run_campaign(ax, label):
    os.makedirs(OUTDIR, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    steps = []


    # Sens + puis sens - : symétrie testée, points doublés.
    amps = [(+t) for t in TAU_LIST] + [(-t) for t in TAU_LIST]

    
    print(f"\n=== RUN {label} : {len(amps)} échelons ===")

    for k, tau in enumerate(amps):
        if not wait_until_stopped(ax):                                                      # Au cas où l'arrêt n'a pas été atteint
            print("    !! arrêt non atteint, on saute cet échelon")
            continue

        print(f"  [{k+1:2d}/{len(amps)}] tau = {tau:+.3f} N.m ...", end="", flush=True)

        s = log_step(ax, tau, SETTLE_HORIZON)
        s["label"] = label
        s["index"] = k

        steps.append(s)                                                                     # On garde les logs dans une liste 

        v_inf = np.median(s["vel"][-max(5, len(s["vel"])//10):])                            # On garde la vitesse de régime (médiane 10% derniers points)
        fs = actual_rate(s)                                                                 
        print(f" theta_dot_inf={v_inf:+6.2f} rad/s   fs={fs:.0f} Hz")                       
        time.sleep(PAUSE)

    # Sauvegarde des logs dans un fichier npz. 
    outfile = os.path.join(OUTDIR, f"id_{label}_{stamp}.npz")
    np.savez(
        outfile,
        label=label,
        tau_cmd=np.array([s["tau_cmd"] for s in steps]),
        t=np.array([s["t"] for s in steps], dtype=object),
        pos=np.array([s["pos"] for s in steps], dtype=object),
        vel=np.array([s["vel"] for s in steps], dtype=object),
        iq=np.array([s["iq"] for s in steps], dtype=object),
    )
    print(f"\nSauvegardé : {outfile}")
    print("Analyse : python3 analyze_identification.py", outfile)
    return outfile


# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Acquisition identification moteur.")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--probe", action="store_true", help="balayage grossier")
    g.add_argument("--run", metavar="LABEL", help="campagne complète (ex: E1a, E1b)")
    args = p.parse_args()

    odrv, ax = connect()
    try:
        enter_torque_mode(ax)
        if args.probe:
            probe(ax)
        else:
            run_campaign(ax, args.run)
    except KeyboardInterrupt:
        print("\nInterrompu.")
    finally:
        stop(ax)
        print("Moteur en IDLE.")


if __name__ == "__main__":
    main()
