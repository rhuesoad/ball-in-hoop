
import argparse
import json
import math
import os
import time
from datetime import datetime

import cv2
import numpy as np
from scipy.linalg import eig

from common.ball_detection import BallDetectorAA4CC, DEFAULT_COLOR_COEFS, open_camera

CAMERA_SIZE = (820, 616)
CAMERA_FPS = 50
EXPOSURE_US = 10000 
GAIN_CAM = 12.0
FIX_CAMERA_CONTROLS = True
AWB_RED_GAIN = 1.8
AWB_BLUE_GAIN = 1.6
FRAME_DURATION_US = (int(1e6 / CAMERA_FPS), int(1e6 / CAMERA_FPS))

THRESHOLD = 90
DOWNSAMPLE = 8
TRACKING_WIN = 128
BALL_SIZE = (0, 150)

OUTER_RADIUS_M = 0.1035
INNER_RADIUS_M = 0.0445
DEFAULT_DURATION = 30.0
DEFAULT_STATIC_DURATION = 3.0

EXPECTED_OUTER_RADIUS_PX = 280.0
OUTER_RADIUS_TOLERANCE_PX = 45.0
EXPECTED_CENTER_PX = (410.0, 310.0)
CENTER_TOLERANCE_PX = 180.0
HOUGH_DP = 1.2
HOUGH_MIN_DIST = 120
HOUGH_PARAM1 = 120
HOUGH_PARAM2 = 50
HOUGH_MIN_RADIUS = int(EXPECTED_OUTER_RADIUS_PX - OUTER_RADIUS_TOLERANCE_PX)
HOUGH_MAX_RADIUS = int(EXPECTED_OUTER_RADIUS_PX + OUTER_RADIUS_TOLERANCE_PX)


def lock_camera(cam):
    if not FIX_CAMERA_CONTROLS:
        return
    controls = {
        "AeEnable": False,
        "AwbEnable": False,
        "ExposureTime": int(EXPOSURE_US),
        "AnalogueGain": float(GAIN_CAM),
        "ColourGains": (float(AWB_RED_GAIN), float(AWB_BLUE_GAIN)),
        "FrameDurationLimits": FRAME_DURATION_US,
    }
    try:
        cam.set_controls(controls)
    except Exception as exc:
        print("[cam] contrôles non verrouillés:", exc)


def detect_hough_circles(frame):
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    gray = cv2.GaussianBlur(gray, (9, 9), 2.0)
    circles = cv2.HoughCircles(
        gray, cv2.HOUGH_GRADIENT, dp=HOUGH_DP,
        minDist=HOUGH_MIN_DIST, param1=HOUGH_PARAM1,
        param2=HOUGH_PARAM2, minRadius=HOUGH_MIN_RADIUS,
        maxRadius=HOUGH_MAX_RADIUS,
    )
    if circles is None:
        return np.empty((0, 3), dtype=float)
    circles = np.asarray(circles[0], dtype=float)
    return circles[np.argsort(circles[:, 2])[::-1]]


def filter_by_center(circles, center, tolerance):
    if len(circles) == 0:
        return circles
    cx, cy = center
    d = np.hypot(circles[:, 0] - cx, circles[:, 1] - cy)
    return circles[d <= tolerance]


def select_expected_circle(circles, expected_radius, tolerance):
    if len(circles) == 0:
        return None
    d = np.abs(circles[:, 2] - expected_radius)
    accepted = circles[d <= tolerance]
    if len(accepted) == 0:
        return None
    return accepted[int(np.argmin(np.abs(accepted[:, 2] - expected_radius)))]


def draw_hough_result(frame, circles, selected, path):
    vis = frame.copy()
    for cx, cy, radius in circles:
        selected_flag = np.allclose(selected, [cx, cy, radius])
        color = (0, 255, 0) if selected_flag else (255, 120, 0)
        thickness = 5 if selected_flag else 2
        cv2.circle(vis, (int(round(cx)), int(round(cy))), int(round(radius)), color, thickness)
        cv2.circle(vis, (int(round(cx)), int(round(cy))), 5, color, -1)
        if selected_flag:
            cv2.putText(vis, f"SELECTED r={radius:.1f}px",
                        (max(5, int(cx - 100)), max(20, int(cy - radius - 10))),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 0), 2)
    cv2.imwrite(path, vis)


def roi_from_circle(cx, cy, radius, margin=0.05):
    width, height = CAMERA_SIZE
    half = radius * (1.0 + margin)
    x0 = max(0, int(round(cx - half)))
    y0 = max(0, int(round(cy - half)))
    x1 = min(width, int(round(cx + half)))
    y1 = min(height, int(round(cy + half)))
    return x0, y0, x1 - x0, y1 - y0


def taubin_circle_fit(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 6:
        raise ValueError("Il faut au moins 6 points 2D.")
    mean = points.mean(axis=0)
    xy = points - mean
    x, y = xy[:, 0], xy[:, 1]
    z = x * x + y * y
    Z = np.column_stack((z, x, y, np.ones_like(x)))
    M = (Z.T @ Z) / len(points)
    N = np.array([[0., 0., 0., 2.], [0., 1., 0., 0.],
                  [0., 0., 1., 0.], [2., 0., 0., 0.]])
    values, vectors = eig(M, N)
    candidates = []
    for value, vector in zip(values, vectors.T):
        if not np.isfinite(value) or abs(value.imag) > 1e-8:
            continue
        q = np.real(vector)
        a, b, c, d = q
        if abs(a) < 1e-12:
            continue
        r2 = (b*b + c*c - 4*a*d) / (4*a*a)
        if r2 > 0 and np.isfinite(r2):
            candidates.append((abs(float(value.real)), q))
    if not candidates:
        raise RuntimeError("Échec de l'ajustement de Taubin.")
    _, q = min(candidates, key=lambda item: item[0])
    a, b, c, d = q
    cx = -b / (2*a) + mean[0]
    cy = -c / (2*a) + mean[1]
    radius = math.sqrt(max((b*b + c*c - 4*a*d) / (4*a*a), 0.0))
    residuals = np.hypot(points[:, 0] - cx, points[:, 1] - cy) - radius
    return float(cx), float(cy), float(radius), residuals


def robust_taubin_fit(points, iterations=4):
    points = np.asarray(points, dtype=float)
    keep = np.ones(len(points), dtype=bool)
    for _ in range(iterations):
        cx, cy, radius, _ = taubin_circle_fit(points[keep])
        all_errors = np.abs(np.hypot(points[:, 0] - cx, points[:, 1] - cy) - radius)
        med = np.median(all_errors)
        mad = np.median(np.abs(all_errors - med))
        keep = all_errors <= med + 3.5 * max(mad, 0.5)
        if keep.sum() < 6:
            break
    cx, cy, radius, residuals = taubin_circle_fit(points[keep])
    return cx, cy, radius, keep, residuals


def draw_taubin_result(frame, points, cx, cy, radius, keep, path):
    vis = frame.copy()
    for point, valid in zip(points, keep):
        color = (0, 255, 0) if valid else (0, 0, 255)
        cv2.circle(vis, tuple(np.round(point).astype(int)), 2, color, -1)
    cv2.circle(vis, (round(cx), round(cy)), round(radius), (255, 0, 255), 3)
    cv2.circle(vis, (round(cx), round(cy)), 6, (255, 0, 255), -1)
    cv2.putText(vis, f"Taubin ({cx:.1f},{cy:.1f}) r={radius:.1f}px",
                (20, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 0, 255), 2)
    cv2.imwrite(path, vis)


def make_detector():
    return BallDetectorAA4CC(
        color_coefs=DEFAULT_COLOR_COEFS,
        threshold=THRESHOLD,
        downsample=DOWNSAMPLE,
        tracking_window=TRACKING_WIN,
        ball_size=BALL_SIZE,
    )


def collect_ball_trajectory(cam, detector, roi, duration):
    rx, ry, rw, rh = roi
    points, times = [], []
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < duration:
        frame = cam.capture_array()
        t = time.perf_counter() - t0
        local = detector.process_image(frame[ry:ry+rh, rx:rx+rw])
        if local is not None:
            points.append((local[0] + rx, local[1] + ry))
            times.append(t)
    return np.asarray(points, float), np.asarray(times, float)


def collect_static_ball(cam, detector, roi, duration, center):
    rx, ry, rw, rh = roi
    cx, cy = center
    angles = []
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < duration:
        frame = cam.capture_array()
        local = detector.process_image(frame[ry:ry+rh, rx:rx+rw])
        if local is not None:
            u, v = local[0] + rx, local[1] + ry
            angles.append(math.atan2(u - cx, v - cy))
    return np.asarray(angles, float)


def main():
    parser = argparse.ArgumentParser(description="Calibration caméra AA4CC")
    parser.add_argument("--duration", type=float, default=DEFAULT_DURATION)
    parser.add_argument("--static-duration", type=float, default=DEFAULT_STATIC_DURATION)
    parser.add_argument("--outdir", default="data/camera_calibration")
    parser.add_argument("--focal-px", type=float, default=680.0)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    prefix = os.path.join(args.outdir, f"camera_calibration_{stamp}")
    cam = None

    try:
        print("[cam] ouverture de la caméra...")
        cam = open_camera(CAMERA_SIZE, CAMERA_FPS, EXPOSURE_US, GAIN_CAM)
        lock_camera(cam)
        detector = make_detector()
        for _ in range(10):
            cam.capture_array()

        frame = cam.capture_array()
        cv2.imwrite(prefix + "_original.png", frame)
        raw_circles = detect_hough_circles(frame)
        if len(raw_circles) == 0:
            raise RuntimeError("Aucun cercle Hough détecté.")

        candidates = filter_by_center(raw_circles, EXPECTED_CENTER_PX, CENTER_TOLERANCE_PX)
        selected = select_expected_circle(candidates, EXPECTED_OUTER_RADIUS_PX, OUTER_RADIUS_TOLERANCE_PX)
        if selected is None:
            raise RuntimeError(
                "Aucun cercle compatible. "
                f"Centre attendu {EXPECTED_CENTER_PX} +/- {CENTER_TOLERANCE_PX}px, "
                f"rayon {EXPECTED_OUTER_RADIUS_PX} +/- {OUTER_RADIUS_TOLERANCE_PX}px."
            )

        hough_cx, hough_cy, hough_radius = map(float, selected)
        roi = roi_from_circle(hough_cx, hough_cy, hough_radius)
        rx, ry, rw, rh = roi
        draw_hough_result(frame, candidates, selected, prefix + "_hough.png")

        print(f"[Hough] centre=({hough_cx:.2f}, {hough_cy:.2f}) px, rayon={hough_radius:.2f} px")
        print(f"[roi] x[{rx}:{rx+rw}] y[{ry}:{ry+rh}] ({rw} x {rh} px)")
        input("Vérifie *_hough.png puis appuie sur ENTREE...")

        print(f"Pendant {args.duration:.1f} s, déplace la balle sur toute la circonférence.")
        input("Appuie sur ENTREE pour commencer...")
        points, t_points = collect_ball_trajectory(cam, detector, roi, args.duration)
        print(f"[ball] détections : {len(points)}")
        if len(points) < 100:
            raise RuntimeError(f"Pas assez de détections ({len(points)}).")

        cx, cy, radius_px, keep, residuals = robust_taubin_fit(points)
        rms_px = float(np.sqrt(np.mean(residuals[keep] ** 2)))
        draw_taubin_result(frame, points, cx, cy, radius_px, keep, prefix + "_taubin.png")
        print(f"[Taubin] centre=({cx:.3f}, {cy:.3f}) px, rayon={radius_px:.3f} px")
        print(f"[Taubin] points retenus={int(keep.sum())}/{len(points)}, RMS={rms_px:.3f} px")

        scale = OUTER_RADIUS_M / radius_px
        distance = args.focal_px * OUTER_RADIUS_M / radius_px if args.focal_px else None
        print(f"[scale] lambda={scale:.9f} m/px")
        print(f"[camera] d_cam={distance:.4f} m" if distance else "[camera] d_cam inconnue")

        input("Place la balle au fond, attends l'immobilité, puis appuie sur ENTREE...")
        static = collect_static_ball(cam, detector, roi, args.static_duration, (cx, cy))
        if len(static):
            offset = float(np.median(static))
            spread = float(np.percentile(static, 99) - np.percentile(static, 1))
            print(f"[static] offset={math.degrees(offset):+.3f} deg, p1-p99={math.degrees(spread):.3f} deg")
        else:
            offset = spread = float("nan")
            print("[static] aucune détection")

        static_frame = cam.capture_array()
        local = detector.process_image(static_frame[ry:ry+rh, rx:rx+rw])
        static_image = static_frame.copy()
        cv2.circle(static_image, (round(cx), round(cy)), round(radius_px), (255, 0, 255), 3)
        cv2.circle(static_image, (round(cx), round(cy)), 6, (255, 0, 255), -1)
        if local is not None:
            u, v = local[0] + rx, local[1] + ry
            cv2.circle(static_image, (round(u), round(v)), 7, (0, 255, 0), -1)
            cv2.line(static_image, (round(cx), round(cy)), (round(u), round(v)), (0, 255, 0), 2)
        cv2.imwrite(prefix + "_static.png", static_image)

        meta = {
            "timestamp": stamp,
            "camera_size": list(CAMERA_SIZE), "camera_fps": CAMERA_FPS,
            "exposure_us": EXPOSURE_US, "gain": GAIN_CAM, "roi": list(roi),
            "aa4cc": {"threshold": THRESHOLD, "downsample": DOWNSAMPLE,
                      "tracking_window": TRACKING_WIN, "ball_size": list(BALL_SIZE)},
            "hough": {"all_circles": raw_circles.tolist(),
                      "filtered_circles": candidates.tolist(),
                      "selected": [hough_cx, hough_cy, hough_radius]},
            "taubin": {"center_px": [cx, cy], "radius_px": radius_px,
                       "n_points": len(points), "n_inliers": int(keep.sum()), "rms_px": rms_px},
            "geometry": {"outer_radius_m": OUTER_RADIUS_M, "lambda_m_per_px": scale,
                         "focal_px": args.focal_px, "camera_distance_m": distance},
            "static": {"offset_deg": math.degrees(offset) if np.isfinite(offset) else None,
                       "p1_p99_deg": math.degrees(spread) if np.isfinite(spread) else None},
        }
        np.savez_compressed(prefix + ".npz", meta=json.dumps(meta),
                            hough_circles_raw=raw_circles,
                            hough_circles_filtered=candidates,
                            hough_selected=selected,
                            trajectory_px=points, trajectory_t=t_points,
                            trajectory_inliers=keep, static_psi=static)
        with open(prefix + ".json", "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)
        print("[out]", prefix + ".npz")
        print(f"HOOP_CENTRE_PX = ({cx:.3f}, {cy:.3f})")
        print(f"HOOP_RADIUS_PX = {radius_px:.3f}")
        print(f"LAMBDA_M_PER_PX = {scale:.9f}")

    finally:
        if cam is not None:
            try:
                cam.stop()
            except Exception:
                pass
            try:
                cam.close()
            except Exception:
                pass


if __name__ == "__main__":
    main()
