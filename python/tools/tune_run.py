#!/usr/bin/env python3
"""
Caracterisation de la boucle de vitesse ODrive par placement de poles.

Objectif : etablir la plus grande pulsation propre wn de la boucle interne
qui reste propre (pas de saturation, dispersion faible, pas d'ondulation
residuelle), afin de justifier la separation d'echelles vis-a-vis de la
dynamique de la balle (w_ball ~ 8.4 rad/s).

Synthese :
    Kp[Nm/(rad/s)]   = 2 * zeta * wn * J
    Ki[Nm/(rad/s)/s] = wn^2 * J
puis conversion en unites ODrive (vitesse en tr/s) par un facteur 2*pi.

Usage :
    python3 tune_run.py --check
    python3 tune_run.py --calib
    python3 tune_run.py --wn 10 --zeta 1.5 --J 2.5e-3 --amp 0.1 -n 5
    python3 tune_run.py --kp 0.47 --ki 1.57 -n 5      (gains absolus directs)
"""

import argparse
import datetime
import pathlib
import threading
import time

import numpy as np
import odrive
from odrive.enums import AxisState, ControlMode, InputMode, ProcedureResult
from odrive.utils import (
    TimestampFmt,
    dump_errors,
    high_rate_capture_start,
)

# Pulsation propre de la balle sur le cerceau exterieur [rad/s].
# w_ball = sqrt(c_o / a_o) avec les coefficients de Gurtner & Zemanek (2017),
# eq. 8, evalues sur la geometrie du banc. C'est la frequence a laquelle la
# separation d'echelles doit etre verifiee.
W_BALL = 8.4

# Candidats de noms de proprietes : le firmware 0.6.x expose la vitesse
# tantot comme raccourci, tantot via le mapper. On resout au demarrage.
VEL_CANDIDATES = ["axis0.vel_estimate", "axis0.pos_vel_mapper.vel"]
FIXED_PROPS = [
    "axis0.controller.input_vel",
    "axis0.motor.foc.Iq_setpoint",
    "axis0.motor.foc.Iq_measured",
]


# ----------------------------------------------------------------------
# Connexion et configuration
# ----------------------------------------------------------------------

def connect():
    odrv = odrive.find_any(timeout=15)
    print(
        f"ODrive {odrv.serial_number:x} | "
        f"firmware {odrv.fw_version_major}."
        f"{odrv.fw_version_minor}."
        f"{odrv.fw_version_revision} | "
        f"bus {odrv.vbus_voltage:.1f} V"
    )
    return odrv


def _getattr_path(root, path):
    obj = root
    for part in path.split("."):
        obj = getattr(obj, part)
    return obj


def resolve_props(odrv):
    """Verifie quel nom de propriete de vitesse existe sur ce firmware."""
    vel_prop = None
    for cand in VEL_CANDIDATES:
        try:
            _getattr_path(odrv, cand)
            vel_prop = cand
            break
        except AttributeError:
            continue

    if vel_prop is None:
        raise RuntimeError(
            f"Aucune propriete de vitesse trouvee parmi {VEL_CANDIDATES}. "
            "Lancer odrivetool et inspecter axis0."
        )

    for prop in FIXED_PROPS:
        _getattr_path(odrv, prop)

    print(f"Signal de vitesse   : {vel_prop}")
    return [vel_prop] + FIXED_PROPS


def show_config(odrv):
    ax = odrv.axis0
    m = ax.config.motor
    c = ax.controller.config

    print("\n--- Configuration ---")
    print(f"vel_gain              : {c.vel_gain:.6g} Nm/(tr/s)")
    print(f"vel_integrator_gain   : {c.vel_integrator_gain:.6g}")
    print(f"encoder_bandwidth     : {ax.config.encoder_bandwidth:.1f} rad/s")
    print(f"current_soft_max      : {m.current_soft_max:.2f} A")
    print(f"current_hard_max      : {m.current_hard_max:.2f} A")
    print(f"torque_constant       : {m.torque_constant:.6g} Nm/A")
    print(f"vel_limit             : {c.vel_limit:.3f} tr/s")
    print(f"current_control_bw    : {m.current_control_bandwidth:.1f} rad/s")
    print(f"control_mode          : {c.control_mode}")
    print(f"input_mode            : {c.input_mode}")

    print("\n--- Erreurs ---")
    dump_errors(odrv)


# ----------------------------------------------------------------------
# Gains
# ----------------------------------------------------------------------

def compute_gains(wn, zeta, J):
    """Placement de poles -> gains en unites ODrive (vitesse en tr/s)."""
    kp_si = 2.0 * zeta * wn * J          # Nm/(rad/s)
    ki_si = wn ** 2 * J                  # Nm/(rad/s)/s
    return kp_si * 2 * np.pi, ki_si * 2 * np.pi


def apply_gains(odrv, kp, ki):
    """Applique des gains ABSOLUS. Retourne les anciennes valeurs."""
    c = odrv.axis0.controller.config
    old_kp, old_ki = c.vel_gain, c.vel_integrator_gain

    c.vel_gain = float(kp)
    c.vel_integrator_gain = float(ki)

    print("\n--- Gains appliques ---")
    print(f"vel_gain            : {old_kp:.6g} -> {c.vel_gain:.6g}")
    print(f"vel_integrator_gain : {old_ki:.6g} -> {c.vel_integrator_gain:.6g}")

    return old_kp, old_ki


# ----------------------------------------------------------------------
# Calibration et armement
# ----------------------------------------------------------------------

def wait_procedure(odrv, timeout=30.0):
    """
    Attend la fin de la procedure en cours (retour a IDLE).

    Une demande de changement d'etat rend la main immediatement : sans cette
    attente on lit procedure_result == BUSY et on conclut a tort a un echec.
    Pire, relancer une procedure pendant qu'une autre tourne l'annule
    (CANCELLED).
    """
    ax = odrv.axis0
    t0 = time.time()

    while time.time() - t0 < timeout:
        if ax.current_state == AxisState.IDLE:
            time.sleep(0.2)
            return ax.procedure_result
        time.sleep(0.1)

    raise RuntimeError(f"Timeout : procedure toujours en cours apres {timeout:.0f} s")


def calibrate(odrv):
    """
    Recalibre l'offset encodeur.

    L'AMT102-V n'a pas d'index : l'offset est perdu a chaque coupure
    d'alimentation de l'ODrive et doit etre refait avant toute boucle fermee.
    On ne chaine PAS automatiquement vers FULL_CALIBRATION_SEQUENCE : celle-ci
    reecrit aussi resistance et inductance du moteur, ce qui doit rester une
    decision explicite.
    """
    ax = odrv.axis0
    odrv.clear_errors()
    time.sleep(0.2)

    print("Calibration offset encodeur (le cerceau va bouger)...")
    ax.requested_state = AxisState.ENCODER_OFFSET_CALIBRATION
    time.sleep(0.5)
    result = wait_procedure(odrv)

    if result != ProcedureResult.SUCCESS:
        dump_errors(odrv)
        raise RuntimeError(
            f"Calibration offset echouee : {result}\n"
            "Si le moteur n'a jamais ete calibre sur cette configuration, "
            "lancer UNE FOIS a la main dans odrivetool :\n"
            "  odrv0.axis0.requested_state = AxisState.FULL_CALIBRATION_SEQUENCE\n"
            "  (attendre le retour a IDLE, puis)\n"
            "  odrv0.save_configuration()"
        )

    print("Calibration OK")


def arm(odrv):
    ax = odrv.axis0

    odrv.clear_errors()
    time.sleep(0.2)

    needs_cal = ax.procedure_result == ProcedureResult.NOT_CALIBRATED
    if not needs_cal:
        try:
            needs_cal = not ax.is_homed
        except AttributeError:
            pass

    if needs_cal:
        calibrate(odrv)

    ax.controller.config.control_mode = ControlMode.VELOCITY_CONTROL
    ax.controller.config.input_mode = InputMode.PASSTHROUGH
    ax.controller.input_vel = 0.0

    ax.requested_state = AxisState.CLOSED_LOOP_CONTROL
    time.sleep(0.5)

    if ax.current_state != AxisState.CLOSED_LOOP_CONTROL:
        dump_errors(odrv)
        raise RuntimeError(f"Armement echoue : etat = {ax.current_state}")


def disarm(odrv):
    ax = odrv.axis0
    try:
        ax.controller.input_vel = 0.0
        time.sleep(0.3)
        ax.requested_state = AxisState.IDLE
    except Exception:
        pass


# ----------------------------------------------------------------------
# Capture
# ----------------------------------------------------------------------

def capture_step(odrv, amplitude, props, delay=0.05, min_pre=10):
    """
    Capture un echelon de vitesse.

    Retourne (t, signals) avec t=0 a l'instant de l'echelon, ou None si la
    capture est invalide. Le rejet min_pre est essentiel : si l'echelon
    survient avant le debut effectif de la capture, l'origine des temps est
    fausse et le t90 mesure n'a aucun sens.
    """
    ax = odrv.axis0
    cap = high_rate_capture_start(odrv, props)

    timer = threading.Timer(
        delay, lambda: setattr(ax.controller, "input_vel", amplitude)
    )
    timer.start()

    try:
        rec = cap.trigger_and_download_sync(
            trigger_point=0.0,
            return_as=np.recarray,
            t_fmt=TimestampFmt.NANOSECONDS,
        )
    finally:
        timer.cancel()
        ax.controller.input_vel = 0.0
        close = getattr(cap, "close", None) or getattr(cap, "stop", None)
        if close is not None:
            try:
                close()
            except Exception:
                pass

    names = rec.dtype.names
    t = np.asarray(rec[names[0]], dtype=float) * 1e-9

    signals = {}
    for name in names[1:]:
        signals[name.split(".")[-1]] = np.asarray(rec[name], dtype=float)

    key_vel = "vel" if "vel" in signals else "vel_estimate"
    signals["velocity_raw"] = signals[key_vel]

    sp = signals["input_vel"]
    indices = np.where(np.abs(sp) >= 0.5 * abs(amplitude))[0]

    if len(indices) == 0:
        print("  REJET : echelon absent de la capture")
        return None

    i0 = indices[0]
    if i0 < min_pre:
        print(
            f"  REJET : echelon a l'indice {i0} (< {min_pre}). "
            "Capture demarree trop tard, origine des temps fausse."
        )
        return None

    return t - t[i0], signals


# ----------------------------------------------------------------------
# Metriques
# ----------------------------------------------------------------------

def metric(t, signals, amplitude, i_limit):
    v = signals["velocity_raw"] / amplitude
    iq = signals["Iq_measured"]

    dt = float(np.median(np.diff(t)))
    post = t >= 0
    tv, vv = t[post], v[post]

    def first_crossing(level):
        idx = np.where(vv >= level)[0]
        return np.nan if len(idx) == 0 else tv[idx[0]]

    def settling(tol=0.05):
        """Dernier instant ou |v-1| sort de la tolerance."""
        out = np.where(np.abs(vv - 1.0) > tol)[0]
        return np.nan if len(out) == 0 else tv[out[-1]]

    # Ondulation residuelle sur le dernier tiers : distingue un transitoire
    # amorti d'un cycle limite entretenu.
    tail = tv >= (tv[0] + 0.66 * (tv[-1] - tv[0]))
    ripple = float(np.std(vv[tail])) * 100.0

    iq_max = float(np.max(np.abs(iq)))
    saturated = iq_max >= 0.95 * i_limit

    print(f"  dt        : {1e6 * dt:.0f} us  ({len(t)} pts)")
    print(f"  t63       : {1e3 * first_crossing(0.632):.1f} ms")
    print(f"  t90       : {1e3 * first_crossing(0.900):.1f} ms")
    print(f"  t_5%      : {1e3 * settling():.1f} ms")
    print(f"  overshoot : {100.0 * (np.max(vv) - 1.0):.1f} %")
    print(f"  ondul.    : {ripple:.2f} %  (dernier tiers)")
    print(
        f"  Iq max    : {iq_max:.2f} A"
        + ("   *** SATURE ***" if saturated else "")
    )

    return {
        "dt": dt,
        "t63": first_crossing(0.632),
        "t90": first_crossing(0.900),
        "t_settle": settling(),
        "overshoot": 100.0 * (np.max(vv) - 1.0),
        "ripple": ripple,
        "iq_max": iq_max,
        "saturated": saturated,
    }


def separation_report(wn, zeta):
    """
    Dephasage theorique apporte par la boucle interne a w_ball.

    Boucle fermee PI + inertie pure :
        T(s) = (2*zeta*wn*s + wn^2) / (s^2 + 2*zeta*wn*s + wn^2)
    Le zero du PI apporte de l'avance de phase : un modele du premier ordre
    surestime largement le retard. C'est le chiffre qui justifie (ou non)
    de negliger theta_dot dans l'etat du LQR externe.
    """
    s = 1j * W_BALL
    T = (2 * zeta * wn * s + wn ** 2) / (s ** 2 + 2 * zeta * wn * s + wn ** 2)
    phase = float(np.degrees(np.angle(T)))
    gain = float(np.abs(T))

    print(f"\n--- Separation d'echelles a w_ball = {W_BALL:.1f} rad/s ---")
    print(f"gain boucle interne   : {gain:.3f}  ({20 * np.log10(gain):+.1f} dB)")
    print(f"dephasage             : {phase:+.1f} deg")
    if phase > -10.0:
        print("=> negligeable : theta_dot peut etre omis de l'etat du LQR.")
    elif phase > -25.0:
        print("=> marginal : documenter le cout en marge de phase.")
    else:
        print("=> non negligeable : augmenter l'etat du LQR de theta_dot.")

    return {"gain": gain, "phase_deg": phase}


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--check", action="store_true", help="affiche la config et sort")
    p.add_argument("--calib", action="store_true", help="calibre l'offset encodeur et sort")
    p.add_argument("--wn", type=float, default=None, help="pulsation propre [rad/s]")
    p.add_argument("--zeta", type=float, default=1.5, help="amortissement")
    p.add_argument("--J", type=float, default=2.5e-3, help="inertie totale [kg.m2]")
    p.add_argument("--kp", type=float, default=None, help="vel_gain absolu (prioritaire)")
    p.add_argument("--ki", type=float, default=None, help="vel_integrator_gain absolu")
    p.add_argument("--amp", type=float, default=0.1, help="amplitude echelon [tr/s]")
    p.add_argument("-n", type=int, default=5, help="nombre d'essais")
    p.add_argument("--restore", action="store_true", help="restaure les gains en sortant")
    args = p.parse_args()

    odrv = connect()
    props = resolve_props(odrv)

    if args.check:
        show_config(odrv)
        return

    if args.calib:
        calibrate(odrv)
        return

    # --- Gains : absolus prioritaires, sinon placement de poles ---
    if args.kp is not None:
        kp = args.kp
        ki = args.ki if args.ki is not None else 0.0
        wn_eff, zeta_eff = np.nan, np.nan
        print(f"\nGains absolus imposes : kp={kp:.6g}, ki={ki:.6g}")
    elif args.wn is not None:
        kp, ki = compute_gains(args.wn, args.zeta, args.J)
        wn_eff, zeta_eff = args.wn, args.zeta
        print(
            f"\nPlacement de poles : wn={args.wn:.2f} rad/s, "
            f"zeta={args.zeta:.2f}, J={args.J:.4g} kg.m2"
        )
    else:
        p.error("Fournir --wn (placement de poles) ou --kp (gain absolu).")

    old_kp, old_ki = apply_gains(odrv, kp, ki)
    show_config(odrv)

    if not np.isnan(wn_eff):
        sep = separation_report(wn_eff, zeta_eff)
    else:
        sep = {"gain": np.nan, "phase_deg": np.nan}

    i_limit = odrv.axis0.config.motor.current_soft_max

    print(f"\nAmplitude : {args.amp:.3f} tr/s ({args.amp * 2 * np.pi:.2f} rad/s)")
    input("Zone degagee, balle retiree. Entree pour commencer... ")

    runs, metrics = [], []

    try:
        arm(odrv)

        for i in range(args.n):
            print(f"\nEssai {i + 1}/{args.n}")
            result = capture_step(odrv, args.amp, props)

            if result is None:
                continue

            t, signals = result
            runs.append((t, signals))
            metrics.append(metric(t, signals, args.amp, i_limit))
            time.sleep(1.0)

        disarm(odrv)

    except KeyboardInterrupt:
        print("\nInterruption utilisateur")
        disarm(odrv)
    except Exception:
        disarm(odrv)
        dump_errors(odrv)
        raise
    finally:
        if args.restore:
            apply_gains(odrv, old_kp, old_ki)

    if not runs:
        raise RuntimeError("Aucun essai valide. Verifier les rejets ci-dessus.")

    # --- Synthese inter-essais ---
    def agg(key):
        vals = np.array([m[key] for m in metrics], dtype=float)
        return np.nanmean(vals), np.nanstd(vals)

    print(f"\n=== Synthese sur {len(runs)} essais valides ===")
    for key, unit, scale in [
        ("t90", "ms", 1e3),
        ("t_settle", "ms", 1e3),
        ("overshoot", "%", 1.0),
        ("ripple", "%", 1.0),
        ("iq_max", "A", 1.0),
    ]:
        mean, std = agg(key)
        print(f"{key:10s}: {scale * mean:7.2f} +/- {scale * std:5.2f} {unit}")

    n_sat = sum(m["saturated"] for m in metrics)
    if n_sat:
        print(f"\nATTENTION : {n_sat}/{len(metrics)} essais saturent en courant.")
        print("Le modele lineaire ne s'applique pas. Reduire --amp.")

    # --- Sauvegarde ---
    outdir = pathlib.Path("data/tuning")
    outdir.mkdir(parents=True, exist_ok=True)

    tag = f"wn{wn_eff:.1f}" if not np.isnan(wn_eff) else "abs"
    fname = outdir / f"run_{tag}_{datetime.datetime.now():%Y%m%d_%H%M%S}.npz"

    payload = {
        "kp": odrv.axis0.controller.config.vel_gain,
        "ki": odrv.axis0.controller.config.vel_integrator_gain,
        "wn": wn_eff,
        "zeta": zeta_eff,
        "J": args.J,
        "amp": args.amp,
        "n_runs": len(runs),
        "w_ball": W_BALL,
        "inner_gain_at_wball": sep["gain"],
        "inner_phase_at_wball_deg": sep["phase_deg"],
        "encoder_bandwidth": odrv.axis0.config.encoder_bandwidth,
        "current_soft_max": i_limit,
        "torque_constant": odrv.axis0.config.motor.torque_constant,
        "current_control_bandwidth": odrv.axis0.config.motor.current_control_bandwidth,
        "vbus": odrv.vbus_voltage,
    }

    for key in ["t90", "t_settle", "overshoot", "ripple", "iq_max", "dt"]:
        payload[f"m_{key}"] = np.array([m[key] for m in metrics], dtype=float)

    for i, (t, s) in enumerate(runs):
        payload[f"t{i}"] = t
        payload[f"velocity{i}"] = s["velocity_raw"]
        payload[f"input_vel{i}"] = s["input_vel"]
        payload[f"iq{i}"] = s["Iq_measured"]
        payload[f"iq_set{i}"] = s["Iq_setpoint"]

    np.savez(fname, **payload)
    print(f"\nDonnees : {fname}")
    print(f"Gains precedents : vel_gain={old_kp:.6g}, vel_integrator_gain={old_ki:.6g}")


if __name__ == "__main__":
    main()
