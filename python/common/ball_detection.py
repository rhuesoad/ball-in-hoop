#!/usr/bin/env python3
"""
ball_detection.py -- Portage fidele du detecteur AA4CC + mesure de latence

PIPELINE 
--------------------------------
    1. Combinaison lineaire des canaux -> image scalaire uint8 :
           im = clip(c0*ch0 + c1*ch1 + c2*ch2, 0, 255)
    2. Masque 
    3. Seuillage : cv2.inRange(im, threshold, 255)
    4. Erosion optionnelle (à checker)
    5. cv2.moments -> centroide (m10/m00, m01/m00)
       Validation de l'aire via m00 (sur image binaire 0/255, m00 = aire*255)

processImage(image) : recherche en deux etages
    a. Decimation image[::d, ::d] -> recherche grossiere
       (limites d'aire divisees par d^2)
    b. Fenetre `tracking_window` centree sur le resultat grossier,
       a pleine resolution -> recherche fine
    c. Transformation des coordonnees ROI -> image complete

Usage
-----
    python3 ball_detection.py                      # 200 frames, defauts
    python3 ball_detection.py --n 500 --fps 50
    python3 ball_detection.py --threshold 90 --preview 5
    python3 ball_detection.py --size 480 480 --exposure-us 10000
"""

import argparse
import json
import math
import os
import sys
import time
from datetime import datetime

import numpy as np
import cv2


# =========================================================================== #
#  PARAMETRES DE REFERENCE AA4CC (config.json_sample, detecteur "Red")
# =========================================================================== #

# Coefficients dans l'ordre des canaux du TABLEAU NUMPY, soit (B, G, R)
# avec picamera2/"RGB888". Equivaut a [1, -0.5, -0.5] en RGB chez AA4CC.
DEFAULT_COLOR_COEFS   = (-0.5, -0.5, 1.0)
DEFAULT_THRESHOLD     = 90      # AA4CC "Red"
DEFAULT_DOWNSAMPLE    = 8
DEFAULT_TRACKING_WIN  = 64
DEFAULT_BALL_SIZE     = (0, 150)  # diametre min/max [px] -> aire via pi*(d/2)^2

DEFAULT_RESOLUTION    = (820, 616)
DEFAULT_FPS           = 50
DEFAULT_EXPOSURE_US   = 10000     # 10 ms, valeur AA4CC
DEFAULT_ANALOGUE_GAIN = 12.0


# =========================================================================== #
#  DETECTEUR (portage de detector.ObjectDetector / BallDetector)
# =========================================================================== #

class BallDetectorAA4CC:
    """Portage fidele de detector.BallDetector (AA4CC / raspi-ballpos)."""

    def __init__(self,
                 color_coefs=DEFAULT_COLOR_COEFS,
                 threshold=DEFAULT_THRESHOLD,
                 downsample=DEFAULT_DOWNSAMPLE,
                 tracking_window=DEFAULT_TRACKING_WIN,
                 ball_size=DEFAULT_BALL_SIZE,
                 mask=None,
                 kernel_dwn=None,
                 kernel_roi=None,
                 debug=False):

        self.color_coefs = tuple(float(c) for c in color_coefs)
        self.threshold = int(threshold)
        self.downsample = int(downsample)
        self.tracking_window = int(tracking_window)
        self.debug = bool(debug)

        # BallDetector.__init__ : ball_size (diametres) -> aires,
        # puis ObjectDetector.__init__ multiplie par 255 car, sur une image
        # binaire 0/255, cv2.moments donne m00 = aire_en_pixels * 255.
        if ball_size is not None:
            a_min = (ball_size[0] / 2.0) ** 2 * math.pi
            a_max = (ball_size[1] / 2.0) ** 2 * math.pi
            self.objectlim = (a_min * 255, a_max * 255)
            self.area_px = (a_min, a_max)
        else:
            self.objectlim = None
            self.area_px = None

        self.mask = mask
        self.mask_dwn = (mask[::self.downsample, ::self.downsample]
                         if mask is not None else None)

        self.kernel_dwn = (np.array(kernel_dwn, np.uint8)
                           if kernel_dwn is not None else None)
        self.kernel_roi = (np.array(kernel_roi, np.uint8)
                           if kernel_roi is not None else None)

        self.location = None
        self.images = {}

    @property
    def tracking_window_halfsize(self):
        return self.tracking_window // 2

    # --------------------------------------------------------------- #
    #  findTheObject : coeur de la detection
    # --------------------------------------------------------------- #
    def find_the_object(self, image, object_size_lim=None, mask=None,
                        name=None, kernel=None, store=False):
        """Retourne (cx, cy) en coordonnees LOCALES a `image`, ou None."""

        # Combinaison lineaire des canaux -> image scalaire "de couleur".
        im = np.clip(
            image[:, :, 0] * self.color_coefs[0] +
            image[:, :, 1] * self.color_coefs[1] +
            image[:, :, 2] * self.color_coefs[2],
            0, 255).astype(np.uint8)

        if store:
            self.images[name] = im

        # Masque optionnel. (Source : cv2.bitwise_and(im, self.im, mask=mask)
        # -- `self.im` n'existe pas, bug corrige ici.)
        if mask is not None:
            im = cv2.bitwise_and(im, im, mask=mask)
            if store:
                self.images[name + "_masked"] = im

        # Seuillage
        im_thrs = cv2.inRange(im, self.threshold, 255)

        if kernel is not None:
            im_thrs = cv2.erode(im_thrs, kernel, iterations=1)

        if store:
            self.images[name + "_thrs"] = im_thrs

        M = cv2.moments(im_thrs)
        if M['m00'] > 0 and (object_size_lim is None or
                             object_size_lim[0] < M['m00'] < object_size_lim[1]):
            return M['m10'] / M['m00'], M['m01'] / M['m00']

        if self.debug:
            print("  [debug] aire hors limites (m00={:.0f}, lim={})".format(
                M['m00'], object_size_lim))
        return None

    # --------------------------------------------------------------- #
    #  processImage : recherche en deux etages
    # --------------------------------------------------------------- #
    def process_image(self, image, store=False):
        """Retourne (x, y) en coordonnees de l'image complete, ou None."""
        h, w = image.shape[:2]

        start_y, end_y = 0, h
        start_x, end_x = 0, w

        if store:
            self.images["image"] = image

        # --- Etage grossier : image decimee ---
        if self.downsample > 1:
            image_dwn = image[::self.downsample, ::self.downsample, :]

            if self.objectlim:
                d2 = self.downsample ** 2
                object_lim_dwn = (self.objectlim[0] // d2,
                                  self.objectlim[1] // d2)
            else:
                object_lim_dwn = None

            location_dwn = self.find_the_object(
                image_dwn, object_size_lim=object_lim_dwn,
                mask=self.mask_dwn, name="downsample",
                kernel=self.kernel_dwn, store=store)

            if not location_dwn:
                if self.debug:
                    print("  [debug] objet absent de l'image complete")
                self.location = None
                return None

            center = (int(self.downsample * location_dwn[0]),
                      int(self.downsample * location_dwn[1]))

            halfsize = self.tracking_window_halfsize
            start_y = max(center[1] - halfsize, 0)
            end_y   = min(center[1] + halfsize, h)
            start_x = max(center[0] - halfsize, 0)
            end_x   = min(center[0] + halfsize, w)

            if store:
                self.images["image_dwn"] = image_dwn

        # --- Etage fin : fenetre a pleine resolution ---
        image_roi = image[start_y:end_y, start_x:end_x, :]
        if store:
            self.images["image_roi"] = image_roi

        mask_roi = (self.mask[start_y:end_y, start_x:end_x]
                    if self.mask is not None else None)

        location_in_roi = self.find_the_object(
            image_roi, object_size_lim=self.objectlim, mask=mask_roi,
            name="roi", kernel=self.kernel_roi, store=store)

        if not location_in_roi:
            if self.debug:
                print("  [debug] objet absent de la ROI")
            self.location = None
            return None

        # ROI -> coordonnees image complete
        x = start_x + location_in_roi[0]
        y = start_y + location_in_roi[1]
        self.location = (x, y)
        return self.location


# =========================================================================== #
#  DETECTION AUTOMATIQUE DU HOOP (une seule fois au demarrage)
# =========================================================================== #

def detect_hoop_roi(image, margin=0.05, verbose=True):
    """
    Localise le hoop exterieur et renvoie le rectangle de recadrage.

    Le hoop est fixe dans l'image : cette detection n'est faite qu'une fois,
    au demarrage, et peut donc se permettre d'etre couteuse.
    """
    h_img, w_img = image.shape[:2]
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    blur = cv2.GaussianBlur(gray, (5, 5), 0)

    circles = cv2.HoughCircles(
        blur, cv2.HOUGH_GRADIENT, dp=1, minDist=h_img // 2,
        param1=80, param2=40,
        minRadius=int(h_img * 0.30), maxRadius=int(h_img * 0.55))

    if circles is not None:
        best = max(circles[0], key=lambda c: c[2])
        cx, cy, r = float(best[0]), float(best[1]), float(best[2])
    else:
        cx, cy, r = w_img / 2.0, h_img / 2.0, h_img * 0.42
        if verbose:
            print("[hoop] Hough n'a rien trouve -> estimation par defaut")

    # --- Raffinement du rayon par scan radial ---
    edges = []
    for i in range(72):
        a = 2 * math.pi * i / 72
        dx, dy = math.cos(a), math.sin(a)
        rs = np.arange(max(0.0, r - 80), r + 80, 0.5)
        xs = (cx + rs * dx).astype(int)
        ys = (cy + rs * dy).astype(int)
        m = (xs >= 0) & (xs < w_img) & (ys >= 0) & (ys < h_img)
        xs, ys, rv = xs[m], ys[m], rs[m]
        if len(xs) < 10:
            continue
        prof = gray[ys, xs].astype(float)
        g = np.abs(np.diff(prof))
        if len(g) == 0:
            continue
        r_edge = rv[int(np.argmax(g))]
        edges.append((cx + r_edge * dx, cy + r_edge * dy, r_edge))

    if len(edges) >= 36:
        e = np.array(edges)
        r = float(np.median(e[:, 2]))
        (cx, cy), _ = cv2.minEnclosingCircle(e[:, :2].astype(np.float32))

    # --- Rectangle de recadrage ---
    half = r * (1.0 + margin)
    x0 = max(int(round(cx - half)), 0)
    y0 = max(int(round(cy - half)), 0)
    x1 = min(int(round(cx + half)), w_img)
    y1 = min(int(round(cy + half)), h_img)

    if verbose:
        lam = 0.237 / (2 * r)
        print("[hoop] centre=({:.0f}, {:.0f})  rayon={:.0f} px  "
              "diametre={:.0f} px".format(cx, cy, r, 2 * r))
        print("[hoop] lambda = {:.3e} m/px | balle 25 mm -> {:.0f} px"
              .format(lam, 0.025 / lam))
        print("[hoop] recadrage : x[{}:{}] y[{}:{}]  ({} x {} px, "
              "marge {:.0f} %)".format(x0, x1, y0, y1,
                                       x1 - x0, y1 - y0, margin * 100))

    return (x0, y0, x1 - x0, y1 - y0)


# =========================================================================== #
#  CAMERA
# =========================================================================== #

def open_camera(size, fps, exposure_us, gain):
    from picamera2 import Picamera2

    picam2 = Picamera2()

    frame_us = int(round(1_000_000 / fps))

    video_config = picam2.create_video_configuration(
        main={
            "size": size,
            "format": "RGB888",
        },
        controls={
            "FrameDurationLimits": (frame_us, frame_us),
            "ExposureTime": int(exposure_us),
            "AnalogueGain": float(gain),
            "AeEnable": False,
            "AwbEnable": False,
            "ColourGains": (1.16, 2.27),
        },
        buffer_count=4,
    )

    picam2.configure(video_config)
    picam2.start()

    print(
        "[cam] {}x{} @ {:.1f} fps | expo {} us | gain {:.1f} | "
        "AWB fixe R={:.2f} B={:.2f}".format(
            size[0], size[1], fps, exposure_us, gain, 1.16, 2.27
        )
    )

    return picam2


# =========================================================================== #
#  BOUCLE PRINCIPALE : acquisition + detection + chronometrage
# =========================================================================== #

def ball_detection(n_frames=200,
                   size=DEFAULT_RESOLUTION,
                   fps=DEFAULT_FPS,
                   exposure_us=DEFAULT_EXPOSURE_US,
                   analogue_gain=DEFAULT_ANALOGUE_GAIN,
                   threshold=DEFAULT_THRESHOLD,
                   color_coefs=DEFAULT_COLOR_COEFS,
                   downsample=DEFAULT_DOWNSAMPLE,
                   tracking_window=DEFAULT_TRACKING_WIN,
                   ball_size=DEFAULT_BALL_SIZE,
                   roi=None,
                   roi_margin=0.05,
                   warmup_frames=10,
                   debug=False,
                   verbose=True):
    """
    Acquiert n_frames, detecte la balle, chronometre chaque etape.
    """
    picam2 = open_camera(size, fps, exposure_us, analogue_gain)

    det = BallDetectorAA4CC(
        color_coefs=color_coefs,
        threshold=threshold,
        downsample=downsample,
        tracking_window=tracking_window,
        ball_size=ball_size,
        debug=debug,
    )

    if verbose:
        print("[det] seuil={} | coefs={} | downsample={} | fenetre={} px"
              .format(threshold, color_coefs, downsample, tracking_window))
        if det.area_px:
            print("[det] aire admise : {:.0f} - {:.0f} px^2"
                  .format(det.area_px[0], det.area_px[1]))

    # Echauffement (frames ignorees)
    for _ in range(warmup_frames):
        picam2.capture_array()

    # --- Recadrage sur le hoop, determine une seule fois ---
    if roi is None:
        roi = detect_hoop_roi(picam2.capture_array(),
                              margin=roi_margin, verbose=verbose)
    rx, ry, rw, rh = roi

    t_acq   = np.empty(n_frames)
    t_det   = np.empty(n_frames)
    t_loop  = np.empty(n_frames)
    t_stamp = np.empty(n_frames)
    U       = np.full(n_frames, np.nan)
    V       = np.full(n_frames, np.nan)
    OK      = np.zeros(n_frames, dtype=bool)

    if verbose:
        print("[run] acquisition de {} frames...".format(n_frames))

    t0 = time.perf_counter()

    for i in range(n_frames):
        ta = time.perf_counter()
        frame = picam2.capture_array()
        tb = time.perf_counter()
      
        loc = det.process_image(frame[ry:ry + rh, rx:rx + rw])
        tc = time.perf_counter()

        t_acq[i]   = (tb - ta) * 1e3
        t_det[i]   = (tc - tb) * 1e3
        t_loop[i]  = (tc - ta) * 1e3
        t_stamp[i] = ta - t0

        if loc is not None:
            # ROI -> coordonnees image complete
            U[i], V[i] = loc[0] + rx, loc[1] + ry
            OK[i] = True

    picam2.stop()
    picam2.close()

    return {
        "t": t_stamp, "u": U, "v": V, "ok": OK,
        "t_acq": t_acq, "t_det": t_det, "t_loop": t_loop,
        "params": {
            "resolution": list(size), "fps": fps,
            "exposure_us": exposure_us, "analogue_gain": analogue_gain,
            "threshold": threshold, "color_coefs": list(color_coefs),
            "downsample": downsample, "tracking_window": tracking_window,
            "ball_size": list(ball_size) if ball_size else None,
            "n_frames": n_frames,
        },
    }


# =========================================================================== #
#  RAPPORT
# =========================================================================== #

def report(res):
    """Affiche le resume des latences et de la detection."""
    t_acq, t_det, t_loop = res["t_acq"], res["t_det"], res["t_loop"]
    ok, t = res["ok"], res["t"]

    def line(label, a):
        print("  {:<22s} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}".format(
            label, a.mean(), np.median(a), a.std(),
            np.percentile(a, 95), a.max()))

    print("\n" + "=" * 72)
    print("  LATENCE")
    print("=" * 72)
    print("  {:<22s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}".format(
        "", "Moy", "Med", "Std", "P95", "Max"))
    print("  " + "-" * 68)
    line("Acquisition [ms]", t_acq)
    line("Detection   [ms]", t_det)
    line("Total       [ms]", t_loop)

    # Cadence reellement obtenue, mesuree sur les intervalles entre frames
    dt = np.diff(t) * 1e3
    print("\n  Intervalle entre frames : med {:.2f} ms  p95 {:.2f} ms"
          .format(np.median(dt), np.percentile(dt, 95)))
    print("  Cadence effective       : {:.1f} fps  (p95 : {:.1f} fps)"
          .format(1e3 / np.median(dt), 1e3 / np.percentile(dt, 95)))
    print("  Part acquisition        : {:.0f} % du temps de boucle"
          .format(100 * t_acq.mean() / t_loop.mean()))

    print("\n  Detection : {}/{} ({:.1f} %)".format(
        int(ok.sum()), ok.size, 100.0 * ok.mean()))

    if ok.sum() >= 2:
        u, v = res["u"][ok], res["v"][ok]
        print("  Position u : {:.1f} - {:.1f} px  (etendue {:.0f})".format(
            u.min(), u.max(), u.max() - u.min()))
        print("  Position v : {:.1f} - {:.1f} px  (etendue {:.0f})".format(
            v.min(), v.max(), v.max() - v.min()))
    print()


def save(res, outdir="data/ball_detection"):
    os.makedirs(outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(outdir, "balldet_{}.npz".format(stamp))
    np.savez_compressed(
        path,
        t=res["t"], u=res["u"], v=res["v"], ok=res["ok"],
        t_acq=res["t_acq"], t_det=res["t_det"], t_loop=res["t_loop"],
        params=json.dumps(res["params"]),
    )
    print("[out] {}".format(path))
    return path


def preview(n, size, fps, exposure_us, analogue_gain, threshold,
            color_coefs, downsample, tracking_window, ball_size,
            roi=None, roi_margin=0.05,
            outdir="data/ball_detection/preview"):
    """Sauvegarde n frames annotees + les images intermediaires du pipeline."""
    picam2 = open_camera(size, fps, exposure_us, analogue_gain)
    det = BallDetectorAA4CC(color_coefs=color_coefs, threshold=threshold,
                            downsample=downsample,
                            tracking_window=tracking_window,
                            ball_size=ball_size,
                            debug=True)
    os.makedirs(outdir, exist_ok=True)

    for _ in range(10):
        picam2.capture_array()

    if roi is None:
        roi = detect_hoop_roi(picam2.capture_array(), margin=roi_margin)
    rx, ry, rw, rh = roi

    for i in range(n):
        frame = picam2.capture_array()
        crop = frame[ry:ry + rh, rx:rx + rw]
        loc = det.process_image(crop, store=True)

        vis = frame.copy()
        cv2.rectangle(vis, (rx, ry), (rx + rw, ry + rh), (0, 255, 255), 1)
        if loc is not None:
            x, y = int(round(loc[0] + rx)), int(round(loc[1] + ry))
            cv2.drawMarker(vis, (x, y), (0, 255, 0), cv2.MARKER_CROSS, 40, 2)
            cv2.putText(vis, "u={} v={}".format(x, y), (10, 25),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 0), 2)
        else:
            cv2.putText(vis, "NON DETECTEE", (10, 25),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 0, 255), 2)

        cv2.imwrite(os.path.join(outdir, "prev_{:02d}.png".format(i)), vis)
        # Projection couleur SUR LE RECADRAGE : c'est ce que voit le seuil.
        proj = np.clip(
            crop[:, :, 0] * color_coefs[0] + crop[:, :, 1] * color_coefs[1] +
            crop[:, :, 2] * color_coefs[2], 0, 255).astype(np.uint8)
        cv2.imwrite(os.path.join(outdir, "prev_{:02d}_proj.png".format(i)), proj)
        cv2.imwrite(os.path.join(outdir, "prev_{:02d}_thrs.png".format(i)),
                    cv2.inRange(proj, threshold, 255))

        print("  frame {} : proj max={} p99={:.0f} | seuil={}".format(
            i, proj.max(), np.percentile(proj, 99), threshold))
        time.sleep(0)

    picam2.stop()
    picam2.close()
    print("[preview] images dans {}".format(outdir))


# =========================================================================== #

def main():
    p = argparse.ArgumentParser(
        description="Detection de balle (portage AA4CC) + mesure de latence")
    p.add_argument("--n", type=int, default=200, help="nombre de frames")
    p.add_argument("--size", type=int, nargs=2, default=list(DEFAULT_RESOLUTION),
                   metavar=("W", "H"))
    p.add_argument("--fps", type=float, default=DEFAULT_FPS)
    p.add_argument("--exposure-us", type=int, default=DEFAULT_EXPOSURE_US)
    p.add_argument("--gain", type=float, default=DEFAULT_ANALOGUE_GAIN)
    p.add_argument("--threshold", type=int, default=DEFAULT_THRESHOLD)
    p.add_argument("--downsample", type=int, default=DEFAULT_DOWNSAMPLE)
    p.add_argument("--tracking-window", type=int, default=DEFAULT_TRACKING_WIN)
    p.add_argument("--ball-size", type=float, nargs=2,
                   default=list(DEFAULT_BALL_SIZE), metavar=("DMIN", "DMAX"),
                   help="diametre min/max de la balle en px")
    p.add_argument("--preview", type=int, metavar="N",
                   help="sauvegarde N frames annotees au lieu de mesurer")
    p.add_argument("--roi-margin", type=float, default=0.05,
                   help="marge du recadrage auto sur le hoop, en fraction "
                        "du rayon (defaut 0.05 = 5%%)")
    p.add_argument("--no-crop", action="store_true",
                   help="desactive le recadrage sur le hoop (image entiere)")
    p.add_argument("--outdir", default="data/ball_detection")
    p.add_argument("--debug", action="store_true")
    args = p.parse_args()

    # --no-crop : ROI = image entiere, donc pas de detection de hoop
    roi = (0, 0, args.size[0], args.size[1]) if args.no_crop else None

    if args.preview:
        preview(args.preview, args.size, args.fps, args.exposure_us,
                args.gain, args.threshold, DEFAULT_COLOR_COEFS,
                args.downsample, args.tracking_window, tuple(args.ball_size),
                roi=roi, roi_margin=args.roi_margin,
                outdir=os.path.join(args.outdir, "preview"))
        return

    res = ball_detection(
        n_frames=args.n, size=args.size, fps=args.fps,
        exposure_us=args.exposure_us, analogue_gain=args.gain,
        threshold=args.threshold, downsample=args.downsample,
        tracking_window=args.tracking_window,
        ball_size=tuple(args.ball_size),
        roi=roi, roi_margin=args.roi_margin, debug=args.debug,
    )
    report(res)
    save(res, args.outdir)


if __name__ == "__main__":
    main()
