#!/usr/bin/env python3

import argparse
import json
import math
import os
import time
from datetime import datetime

import cv2
import numpy as np
import odrive
from odrive.enums import *

from common.ball_detection import (
    BallDetectorAA4CC,
    DEFAULT_COLOR_COEFS,
    open_camera,
)


# ============================================================================
# CAMERA ET DETECTION
# ============================================================================

SIZE = (820, 616)
FPS = 50
EXPOSURE_US = 18000
GAIN = 12.0

THRESHOLD = 110
DOWNSAMPLE = 8
TRACKING_WIN = 64
BALL_SIZE = (0, 150)

ROI_MARGIN = 0.05

# Paramètres Hough déterminés à partir de l'essai précédent
EXPECTED_CENTER_PX = (410.0, 310.0)
CENTER_TOLERANCE_PX = 180.0

EXPECTED_OUTER_RADIUS_PX = 280.0
OUTER_RADIUS_TOLERANCE_PX = 45.0

HOUGH_DP = 1.2
HOUGH_MIN_DIST = 120
HOUGH_PARAM1 = 120
HOUGH_PARAM2 = 50
HOUGH_MIN_RADIUS = int(
    EXPECTED_OUTER_RADIUS_PX - OUTER_RADIUS_TOLERANCE_PX
)
HOUGH_MAX_RADIUS = int(
    EXPECTED_OUTER_RADIUS_PX + OUTER_RADIUS_TOLERANCE_PX
)

CAM_LATENCY_S = 0.0193


# ============================================================================
# MOTEUR
# ============================================================================

TORQUE_CONSTANT = 0.031


# ============================================================================
# HOUGH ET ROI
# ============================================================================

def detect_hough_circles(frame):
    """Détecte les cercles dans l'image complète."""
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
    """Conserve les cercles proches du centre attendu."""
    if len(circles) == 0:
        return circles

    cx_expected, cy_expected = expected_center

    distance = np.hypot(
        circles[:, 0] - cx_expected,
        circles[:, 1] - cy_expected,
    )

    return circles[distance <= tolerance]


def select_expected_circle(
    circles,
    expected_radius,
    tolerance,
):
    """Sélectionne le rayon le plus proche du rayon attendu."""
    if len(circles) == 0:
        return None

    distance = np.abs(circles[:, 2] - expected_radius)
    accepted = circles[distance <= tolerance]

    if len(accepted) == 0:
        return None

    index = int(np.argmin(
        np.abs(accepted[:, 2] - expected_radius)
    ))

    return accepted[index]


def roi_from_circle(
    cx,
    cy,
    radius,
    margin=ROI_MARGIN,
):
    """Construit une ROI autour du cercle Hough."""
    image_width, image_height = SIZE

    half_size = radius * (1.0 + margin)

    x0 = max(0, int(round(cx - half_size)))
    y0 = max(0, int(round(cy - half_size)))
    x1 = min(image_width, int(round(cx + half_size)))
    y1 = min(image_height, int(round(cy + half_size)))

    return (
        x0,
        y0,
        x1 - x0,
        y1 - y0,
    )


def draw_hough_result(
    frame,
    circles,
    selected,
    path,
):
    """Sauvegarde une image avec les candidats Hough."""
    image = frame.copy()

    for cx, cy, radius in circles:
        is_selected = np.allclose(
            selected,
            [cx, cy, radius],
        )

        color = (0, 255, 0) if is_selected else (255, 120, 0)
        thickness = 5 if is_selected else 2

        cv2.circle(
            image,
            (int(round(cx)), int(round(cy))),
            int(round(radius)),
            color,
            thickness,
        )

        cv2.circle(
            image,
            (int(round(cx)), int(round(cy))),
            5,
            color,
            -1,
        )

        if is_selected:
            cv2.putText(
                image,
                f"SELECTED r={radius:.1f}px",
                (
                    max(5, int(cx - 100)),
                    max(20, int(cy - radius - 10)),
                ),
                cv2.FONT_HERSHEY_SIMPLEX,
                0.7,
                (0, 255, 0),
                2,
            )

    cv2.imwrite(path, image)


def detect_hoop_hough(frame, output_path=None):
    """
    Détecte le hoop extérieur et retourne :
        centre_x, centre_y, rayon, roi, cercles_bruts
    """
    circles = detect_hough_circles(frame)

    if len(circles) == 0:
        raise RuntimeError(
            "Aucun cercle Hough détecté."
        )

    candidates = filter_by_center(
        circles,
        EXPECTED_CENTER_PX,
        CENTER_TOLERANCE_PX,
    )

    selected = select_expected_circle(
        candidates,
        EXPECTED_OUTER_RADIUS_PX,
        OUTER_RADIUS_TOLERANCE_PX,
    )

    if selected is None:
        print("[Hough] cercles détectés :")

        for i, (cx, cy, radius) in enumerate(circles):
            print(
                f"  {i}: centre=({cx:.1f},{cy:.1f}), "
                f"rayon={radius:.1f}px"
            )

        raise RuntimeError(
            "Aucun cercle compatible avec les contraintes Hough. "
            f"Centre attendu={EXPECTED_CENTER_PX}, "
            f"rayon attendu={EXPECTED_OUTER_RADIUS_PX:.1f}px."
        )

    cx, cy, radius = map(float, selected)
    roi = roi_from_circle(cx, cy, radius)

    if output_path is not None:
        draw_hough_result(
            frame,
            candidates,
            selected,
            output_path,
        )

    return cx, cy, radius, roi, circles


# ============================================================================
# ODRIVE
# ============================================================================

def preflight_check(odrv0):
    """Vérifications avant mouvement."""
    try:
        brake = odrv0.config.brake_resistor0

        if not brake.enable:
            raise RuntimeError(
                "La résistance de freinage est désactivée."
            )

        print(
            "[odrv] frein : "
            f"{brake.resistance:.1f} ohm, active"
        )
    except AttributeError:
        print(
            "[odrv] impossible de lire brake_resistor0 ; "
            "vérifie manuellement la résistance de freinage."
        )

    print(
        "[odrv] Vbus = "
        f"{odrv0.vbus_voltage:.2f} V"
    )


def check_errors(ax, context=""):
    if ax.active_errors != 0 or ax.disarm_reason != 0:
        raise RuntimeError(
            "Erreur ODrive "
            f"{context}: active=0x{ax.active_errors:X}, "
            f"disarm=0x{ax.disarm_reason:X}"
        )


def connect_and_prepare(
    vel_limit_tps,
    accel_limit_tps2,
    calibrate=False,
):
    """
    Connexion ODrive.

    --calibrate :
        réalise la calibration moteur/encodeur une seule fois.

    Sans --calibrate :
        aucune calibration n'est lancée ; l'axe doit déjà avoir une
        référence encodeur valide.
    """
    print("[odrv] connexion...")
    odrv0 = odrive.find_any()
    ax = odrv0.axis0

    print(
        "[odrv] connecté : "
        f"serial={odrv0.serial_number}, "
        f"fw={odrv0.fw_version_major}."
        f"{odrv0.fw_version_minor}."
        f"{odrv0.fw_version_revision}"
    )

    preflight_check(odrv0)

    if calibrate:
        print("[odrv] calibration moteur + encodeur...")
        ax.requested_state = AXIS_STATE_FULL_CALIBRATION_SEQUENCE

        while ax.current_state != AXIS_STATE_IDLE:
            time.sleep(0.1)

        check_errors(ax, "après calibration")

        if ax.procedure_result != PROCEDURE_RESULT_SUCCESS:
            raise RuntimeError(
                "Calibration échouée : "
                f"procedure_result={ax.procedure_result}"
            )

        print("[odrv] calibration terminée.")
    else:
        print("[odrv] aucune recalibration demandée.")

    check_errors(ax, "avant boucle fermée")

    if ax.current_state != AXIS_STATE_CLOSED_LOOP_CONTROL:
        print(
            "[odrv] passage en boucle fermée "
            "sans recalibration..."
        )

        ax.requested_state = AXIS_STATE_CLOSED_LOOP_CONTROL
        time.sleep(0.5)

    if ax.current_state != AXIS_STATE_CLOSED_LOOP_CONTROL:
        raise RuntimeError(
            "Impossible de passer en boucle fermée sans "
            "recalibrer. Utilise --calibrate une seule fois."
        )

    check_errors(ax, "après boucle fermée")

    ax.controller.config.control_mode = (
        CONTROL_MODE_POSITION_CONTROL
    )
    ax.controller.config.input_mode = (
        INPUT_MODE_TRAP_TRAJ
    )

    ax.controller.config.vel_limit = max(
        ax.controller.config.vel_limit,
        1.5 * vel_limit_tps,
    )

    ax.trap_traj.config.vel_limit = vel_limit_tps
    ax.trap_traj.config.accel_limit = accel_limit_tps2
    ax.trap_traj.config.decel_limit = accel_limit_tps2

    print(
        "[odrv] TRAP_TRAJ : "
        f"vel={vel_limit_tps:.2f} turn/s, "
        f"accel={accel_limit_tps2:.2f} turn/s²"
    )

    return odrv0, ax


def wait_until_still(ax, timeout=3.0):
    start = time.monotonic()

    while abs(ax.vel_estimate) > 0.05:
        if time.monotonic() - start > timeout:
            print("[odrv] timeout wait_until_still")
            break

        time.sleep(0.02)


# ============================================================================
# ACQUISITION
# ============================================================================

def run(
    step_deg,
    hoop,
    accel_tps2,
    vel_tps,
    duration_s,
    baseline_s,
    countdown_s,
    trial,
    outdir,
    calibrate=False,
):
    picam2 = None
    ax = None

    try:
        # ------------------------------------------------------------
        # Camera
        # ------------------------------------------------------------
        print("[cam] ouverture...")
        picam2 = open_camera(
            SIZE,
            FPS,
            EXPOSURE_US,
            GAIN,
        )

        det = BallDetectorAA4CC(
            color_coefs=DEFAULT_COLOR_COEFS,
            threshold=THRESHOLD,
            downsample=DOWNSAMPLE,
            tracking_window=TRACKING_WIN,
            ball_size=BALL_SIZE,
        )

        for _ in range(10):
            picam2.capture_array()

        # ------------------------------------------------------------
        # Hough : centre et ROI
        # ------------------------------------------------------------
        frame_hough = picam2.capture_array()

        os.makedirs(outdir, exist_ok=True)

        hough_preview = os.path.join(
            outdir,
            "last_hough_preview.png",
        )

        cx, cy, hough_radius, roi, hough_circles = (
            detect_hoop_hough(
                frame_hough,
                hough_preview,
            )
        )

        rx, ry, rw, rh = roi

        print(
            "[e3] centre Hough : "
            f"({cx:.2f}, {cy:.2f}) px"
        )
        print(
            "[e3] rayon Hough : "
            f"{hough_radius:.2f} px"
        )
        print(
            "[e3] ROI : "
            f"x[{rx}:{rx + rw}] "
            f"y[{ry}:{ry + rh}] "
            f"({rw} x {rh} px)"
        )

        # ------------------------------------------------------------
        # ODrive
        # ------------------------------------------------------------
        _, ax = connect_and_prepare(
            vel_tps,
            accel_tps2,
            calibrate=calibrate,
        )

        theta0_turns = ax.pos_estimate
        step_turns = math.radians(step_deg) / (2.0 * math.pi)
        target_turns = theta0_turns + step_turns

        print(
            "[e3] theta0="
            f"{theta0_turns:.4f} turn, "
            f"échelon={step_deg:.1f} deg"
        )

        # ------------------------------------------------------------
        # Compte à rebours
        # ------------------------------------------------------------
        print()

        for k in range(int(countdown_s), 0, -1):
            print(
                "  balle au repos : "
                f"départ dans {k}..."
            )
            time.sleep(1.0)

        print(
            f"  >>> ACQUISITION ({duration_s:.1f} s) <<<"
        )

        # ------------------------------------------------------------
        # Buffers
        # ------------------------------------------------------------
        n_max = int(duration_s * FPS * 1.5) + 100

        t_cam = np.full(n_max, np.nan)
        t_enc = np.full(n_max, np.nan)
        theta = np.full(n_max, np.nan)
        theta_dot = np.full(n_max, np.nan)
        iq = np.full(n_max, np.nan)
        u_data = np.full(n_max, np.nan)
        v_data = np.full(n_max, np.nan)
        psi = np.full(n_max, np.nan)
        ok = np.zeros(n_max, dtype=bool)

        step_time = np.nan
        step_sent = False
        i = 0
        t0 = time.perf_counter()

        try:
            while True:
                frame = picam2.capture_array()

                t_cam_i = (
                    time.perf_counter()
                    - t0
                    - CAM_LATENCY_S
                )

                pos_turns = ax.pos_estimate
                vel_turns_s = ax.vel_estimate
                iq_i = ax.motor.foc.Iq_measured
                t_enc_i = time.perf_counter() - t0

                if (
                    not step_sent
                    and t_enc_i >= baseline_s
                ):
                    ax.controller.input_pos = target_turns
                    step_time = time.perf_counter() - t0
                    step_sent = True

                    print(
                        "[e3] échelon envoyé à "
                        f"t={step_time:.3f} s"
                    )

                if (
                    t_enc_i >= duration_s
                    or i >= n_max
                ):
                    break

                loc = det.process_image(
                    frame[ry:ry + rh, rx:rx + rw]
                )

                t_cam[i] = t_cam_i
                t_enc[i] = t_enc_i
                theta[i] = (
                    pos_turns - theta0_turns
                ) * 2.0 * math.pi
                theta_dot[i] = (
                    vel_turns_s * 2.0 * math.pi
                )
                iq[i] = iq_i

                if loc is not None:
                    u, v = loc[0] + rx, loc[1] + ry
                    u_data[i] = u
                    v_data[i] = v

                    # Centre Hough, et non centre de la ROI
                    psi[i] = math.atan2(
                        u - cx,
                        v - cy,
                    )
                    ok[i] = True

                i += 1

                check_errors(ax, "pendant acquisition")

                if i % FPS == 0:
                    detected = (
                        f"{math.degrees(psi[i - 1]):.1f} deg"
                        if ok[i - 1]
                        else "--"
                    )

                    print(
                        f"\r[e3] t={t_enc_i:5.1f}s "
                        f"theta={math.degrees(theta[i - 1]):6.1f}° "
                        f"psi={detected:>10s} "
                        f"Iq={iq_i:6.3f} A",
                        end="",
                        flush=True,
                    )

            print()

        finally:
            print("[e3] retour à l'origine...")
            ax.controller.input_pos = theta0_turns
            time.sleep(1.5)
            wait_until_still(ax)
            ax.requested_state = AXIS_STATE_IDLE

        # ------------------------------------------------------------
        # Découpage
        # ------------------------------------------------------------
        sl = slice(0, i)

        t_cam = t_cam[sl]
        t_enc = t_enc[sl]
        theta = theta[sl]
        theta_dot = theta_dot[sl]
        iq = iq[sl]
        u_data = u_data[sl]
        v_data = v_data[sl]
        psi = psi[sl]
        ok = ok[sl]

        dt = np.diff(t_enc)
        fps_eff = (
            1.0 / np.median(dt)
            if len(dt)
            else float("nan")
        )

        detection_rate = (
            100.0 * ok.mean()
            if len(ok)
            else 0.0
        )

        print(
            f"[e3] {i} échantillons, "
            f"{fps_eff:.1f} Hz effectifs"
        )
        print(
            f"[e3] détection balle : "
            f"{detection_rate:.1f} %"
        )

        if i:
            theta_final = math.degrees(theta[-1])
            print(
                f"[e3] theta final : "
                f"{theta_final:.1f} deg"
            )

        # ------------------------------------------------------------
        # Sauvegarde
        # ------------------------------------------------------------
        stamp = datetime.now().strftime(
            "%Y%m%d_%H%M%S"
        )

        filename = (
            f"E3_{hoop}_step{step_deg:02.0f}deg_"
            f"a{accel_tps2:02.0f}_t{trial:02d}_"
            f"{stamp}.npz"
        )

        path = os.path.join(outdir, filename)

        metadata = {
            "experiment": "E3",
            "trial": trial,
            "timestamp": stamp,
            "hoop": hoop,
            "step_deg": step_deg,
            "t_step_s": (
                float(step_time)
                if np.isfinite(step_time)
                else None
            ),
            "baseline_s": baseline_s,
            "duration_s": duration_s,
            "calibration_performed": bool(calibrate),
            "trap_traj": {
                "accel_limit_tps2": accel_tps2,
                "vel_limit_tps": vel_tps,
                "accel_limit_rad_s2": (
                    accel_tps2 * 2.0 * math.pi
                ),
            },
            "theta0_turns": float(theta0_turns),
            "cam_latency_s": CAM_LATENCY_S,
            "camera": {
                "size": list(SIZE),
                "fps": FPS,
                "exposure_us": EXPOSURE_US,
                "gain": GAIN,
            },
            "hough": {
                "center_px": [cx, cy],
                "radius_px": hough_radius,
                "expected_center_px": list(
                    EXPECTED_CENTER_PX
                ),
                "expected_radius_px": (
                    EXPECTED_OUTER_RADIUS_PX
                ),
                "all_circles": hough_circles.tolist(),
            },
            "roi": [rx, ry, rw, rh],
            "detector": {
                "threshold": THRESHOLD,
                "color_coefs": list(
                    DEFAULT_COLOR_COEFS
                ),
                "downsample": DOWNSAMPLE,
                "tracking_window": TRACKING_WIN,
                "ball_size": list(BALL_SIZE),
            },
            "sample_rate_hz_effective": float(
                fps_eff
            ),
            "detection_rate_pct": float(
                detection_rate
            ),
            "torque_constant": TORQUE_CONSTANT,
        }

        np.savez_compressed(
            path,
            t_cam=t_cam,
            t_enc=t_enc,
            theta=theta,
            theta_dot=theta_dot,
            iq=iq,
            u=u_data,
            v=v_data,
            psi=psi,
            ok=ok,
            meta=json.dumps(metadata),
        )

        print("[out]", path)
        print("[out] aperçu Hough :", hough_preview)

        return path

    finally:
        if picam2 is not None:
            try:
                picam2.stop()
            except Exception:
                pass

            try:
                picam2.close()
            except Exception:
                pass

        if ax is not None:
            try:
                if ax.current_state != AXIS_STATE_IDLE:
                    ax.requested_state = AXIS_STATE_IDLE
            except Exception:
                pass


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="E3 : réponse à un échelon de position"
    )

    parser.add_argument(
        "--step",
        type=float,
        default=45.0,
        help="amplitude de l'échelon en degrés",
    )

    parser.add_argument(
        "--accel",
        type=float,
        default=10.0,
        help="accélération en turn/s²",
    )

    parser.add_argument(
        "--vel",
        type=float,
        default=2.0,
        help="vitesse maximale en turn/s",
    )

    parser.add_argument(
        "--duration",
        type=float,
        default=8.0,
        help="durée totale en secondes",
    )

    parser.add_argument(
        "--baseline",
        type=float,
        default=0.5,
        help="durée avant l'échelon en secondes",
    )

    parser.add_argument(
        "--countdown",
        type=float,
        default=3.0,
    )

    parser.add_argument(
        "--trial",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--outdir",
        default="data/e3",
    )

    parser.add_argument(
        "--hoop",
        choices=["outer", "inner"],
        default="outer",
    )

    parser.add_argument(
        "--calibrate",
        action="store_true",
        help=(
            "calibre moteur et encodeur une seule fois "
            "au démarrage"
        ),
    )

    args = parser.parse_args()

    slip_limit = (
        21.0
        if args.hoop == "outer"
        else 49.0
    )

    if args.accel > slip_limit:
        print(
            "[!] accélération au-dessus de la limite "
            "indicative de roulement sans glissement : "
            f"{args.accel:.1f} > {slip_limit:.1f} turn/s²"
        )

    run(
        step_deg=args.step,
        hoop=args.hoop,
        accel_tps2=args.accel,
        vel_tps=args.vel,
        duration_s=args.duration,
        baseline_s=args.baseline,
        countdown_s=args.countdown,
        trial=args.trial,
        outdir=args.outdir,
        calibrate=args.calibrate,
    )


if __name__ == "__main__":
    main()