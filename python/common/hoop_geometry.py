#!/usr/bin/env python3
"""
hoop_geometry_v3.py -- mesure du centre et du rayon d'orbite de la balle.

POURQUOI CE FICHIER EXISTE
==========================
psi est calcule par control_loop comme

    psi = atan2(u_px - cx, v_px - cy)

ou (cx, cy) = HOOP_CENTRE_PX, une constante figee dans t1_config_v3.py. Tant
que la balle reste au fond (T1/T2), une erreur sur ce centre est absorbee par
l'offset mesure a l'armement : elle est constante et se retranche. Sur un tour
complet elle ne se retranche plus. Un decalage (dx, dy) du centre produit sur
psi une erreur

    dpsi(psi) ~= ( dx cos(psi) - dy sin(psi) ) / R_px

soit une modulation SINUSOIDALE d'amplitude hypot(dx, dy) / R_px, nulle en deux
points seulement et de signe oppose de part et d'autre. Pour 25 px d'ecart sur
une orbite de 208 px cela fait +/- 6.9 deg, et psi_dot en herite d'une
modulation du meme ordre relatif. L'EKF, lui, croit mesurer un angle propre :
aucun reglage de Q ou R ne rattrape un biais qui suit la trajectoire.

Les journaux de la campagne du 11/08 montrent le centre Hough varier de
(385, 275) a (419, 324) d'un essai a l'autre, contre (384.6, 276.6) configure.
Hough n'est donc pas une source de verite exploitable ici : il est refait a
chaque lancement, sur une image, avec un eclairage et une position de camera
qui bougent.

CE QUE FAIT CE SCRIPT
=====================
Il mesure la geometrie UNE FOIS, sur la trajectoire reelle de la balle, et non
sur le cercle apparent du cerceau :

  1. l'operateur fait rouler la balle A LA MAIN, lentement, sur un tour complet
     (plusieurs tours valent mieux) ;
  2. le detecteur journalise le centre de la balle image par image ;
  3. un cercle est ajuste sur le nuage de points au sens des moindres carres
     GEOMETRIQUES (distance orthogonale), ce qui donne d'un seul coup
        - le centre en pixels          -> HOOP_CENTRE_PX
        - le rayon d'ORBITE en pixels  -> HOOP_RADIUS_PX
        - l'echelle lambda             -> R_EFF_M / R_px
        - le residu, donc la confiance qu'on peut leur accorder.

C'est la seule mesure qui repond a la question posee : HOOP_RADIUS_PX doit
decrire le cercle decrit par le CENTRE de la balle, pas le cerceau. Les deux
different du rayon de roulement, et Hough ne mesure ni l'un ni l'autre.

METHODE D'AJUSTEMENT
====================
Deux etapes, la seconde parce que la premiere est biaisee :

  Kasa (algebrique). On resout au sens des moindres carres lineaires
      2 a x_i + 2 b y_i + c = x_i^2 + y_i^2
  d'ou centre (a, b) et R = sqrt(c + a^2 + b^2). Solution exacte et sans
  reglage, mais elle minimise un residu ALGEBRIQUE, pas la distance au cercle :
  l'estimateur est biaise, d'autant plus que l'arc couvert est court
  [Kasa, IEEE Trans. Instrum. Meas. 25(1):8-14, 1976 ; Chernov & Lesort,
   J. Math. Imaging Vis. 23(3):239-252, 2005, sect. 2].

  Gauss-Newton (geometrique). On raffine ensuite en minimisant le vrai cout
      sum_i ( sqrt((x_i-a)^2 + (y_i-b)^2) - R )^2
  dont le jacobien est analytique :
      dr/da = -(x-a)/d,  dr/db = -(y-b)/d,  dr/dR = -1.
  C'est l'estimateur des moindres carres orthogonaux, non biaise a l'ordre
  dominant [Chernov & Lesort, ibid., sect. 3]. Sur un tour complet les deux
  coincident a une fraction de pixel ; l'ecart entre les deux est affiche,
  car un ecart notable signale une couverture angulaire insuffisante.

La couverture angulaire est verifiee explicitement : un cercle ajuste sur un
arc de moins de ~180 deg est mal conditionne, quelle que soit la methode.

SORTIE
======
  - les valeurs a recopier dans la configuration, affichees ;
  - un fichier JSON de calibration, relu par check_geometry() ci-dessous et
    appele par t3_run_v3.py : une configuration qui ne correspond plus a la
    derniere mesure geometrique refuse de partir, au lieu de produire un essai
    ininterpretable.

Usage
-----
    T1_HOOP=outer python3 hoop_geometry_v3.py                 # 30 s
    T1_HOOP=outer python3 hoop_geometry_v3.py --duration 45
    T1_HOOP=inner python3 hoop_geometry_v3.py --exposure 18000

Aucun acces a l'ODrive : le cerceau ne bouge pas, c'est l'operateur qui roule
la balle. Le script est donc sans risque materiel.
"""

import argparse
import json
import math
import os
import time

import numpy as np

from common.ball_detection import (BallDetectorAA4CC, open_camera,
                            DEFAULT_COLOR_COEFS)


CALIB_DIR = os.path.dirname(os.path.abspath(__file__))


def calib_path(hoop):
    return os.path.join(CALIB_DIR, "hoop_geometry_{}.json".format(hoop))


# =========================================================================== #
#  Ajustement de cercle
# =========================================================================== #

def fit_circle_kasa(x, y):
    """Ajustement algebrique. Retourne (cx, cy, R).

    Systeme lineaire  M @ [2a, 2b, c] = x^2 + y^2, resolu par lstsq : pas
    d'inversion explicite, donc robuste a un mauvais conditionnement.
    """
    M = np.column_stack([x, y, np.ones(x.size)])
    rhs = x ** 2 + y ** 2
    sol, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    a, b = 0.5 * sol[0], 0.5 * sol[1]
    R = math.sqrt(max(sol[2] + a ** 2 + b ** 2, 0.0))
    return a, b, R


def fit_circle_geometric(x, y, p0, n_iter=50, tol=1e-9):
    """Moindres carres orthogonaux par Gauss-Newton. Retourne (cx, cy, R).

    Le residu r_i = d_i - R est la distance SIGNEE au cercle. Le jacobien est
    exact ; aucune approximation numerique, aucun pas a regler.
    """
    a, b, R = p0
    for _ in range(n_iter):
        dx, dy = x - a, y - b
        d = np.hypot(dx, dy)
        d = np.where(d < 1e-12, 1e-12, d)
        r = d - R
        J = np.column_stack([-dx / d, -dy / d, -np.ones(d.size)])
        step, *_ = np.linalg.lstsq(J, -r, rcond=None)
        a, b, R = a + step[0], b + step[1], R + step[2]
        if np.max(np.abs(step)) < tol:
            break
    return a, b, R

def fit_circle_robust(
    x,
    y,
    k_sigma=3.0,
    n_iter=8,
    min_keep=200,
    min_sigma=1.0,
):
    """
    Ajustement de cercle avec élagage itératif des points aberrants.

    Retourne :
        cx, cy, R,
        mask_inliers,
        historique
    """
    if x.size < min_keep:
        raise ValueError(
            "pas assez de points pour l'ajustement robuste : {}".format(
                x.size
            )
        )

    mask = np.ones(x.size, dtype=bool)
    historique = []

    for iteration in range(n_iter):
        xi = x[mask]
        yi = y[mask]

        a0, b0, R0 = fit_circle_kasa(xi, yi)
        cx, cy, R = fit_circle_geometric(
            xi, yi, (a0, b0, R0)
        )

        residuals = circle_residuals(x, y, cx, cy, R)

        # Échelle robuste : MAD, non la moyenne et l'écart-type.
        # L'écart-type serait fortement gonflé par les faux points.
        med = float(np.median(residuals))
        mad = float(np.median(np.abs(residuals - med)))
        robust_sigma = max(1.4826 * mad, min_sigma)

        threshold = max(k_sigma * robust_sigma, 2.0)
        new_mask = np.abs(residuals - med) <= threshold

        # Ne jamais supprimer presque tous les points.
        if np.count_nonzero(new_mask) < min_keep:
            order = np.argsort(np.abs(residuals - med))
            new_mask = np.zeros_like(mask)
            new_mask[order[:min_keep]] = True

        n_rejected = int(np.count_nonzero(~new_mask))

        historique.append({
            "iteration": iteration + 1,
            "n_keep": int(np.count_nonzero(new_mask)),
            "n_rejected": n_rejected,
            "sigma_robuste_px": robust_sigma,
            "seuil_px": threshold,
            "rms_tous_points_px": float(
                np.sqrt(np.mean(residuals ** 2))
            ),
        })

        if np.array_equal(new_mask, mask):
            mask = new_mask
            break

        mask = new_mask

    # Ajustement final uniquement sur les inliers.
    xi = x[mask]
    yi = y[mask]

    a0, b0, R0 = fit_circle_kasa(xi, yi)
    cx, cy, R = fit_circle_geometric(
        xi, yi, (a0, b0, R0)
    )

    return cx, cy, R, mask, historique

def circle_residuals(x, y, cx, cy, R):
    return np.hypot(x - cx, y - cy) - R


def angular_coverage(x, y, cx, cy, n_sectors=24):
    """Fraction des secteurs angulaires visites, et plus grand trou en degres.

    Un cercle ajuste sur un arc court est mal conditionne : centre et rayon
    deviennent fortement correles et le residu cesse d'etre informatif.
    """
    ang = np.arctan2(y - cy, x - cx)
    idx = np.floor((ang + math.pi) / (2 * math.pi) * n_sectors).astype(int)
    idx = np.clip(idx, 0, n_sectors - 1)
    seen = np.zeros(n_sectors, dtype=bool)
    seen[idx] = True

    # Plus long trou circulaire de secteurs consecutifs non visites.
    if seen.all():
        return 1.0, 0.0
    doubled = np.concatenate([seen, seen])
    longest = run = 0
    for s in doubled:
        run = 0 if s else run + 1
        longest = max(longest, run)
    longest = min(longest, n_sectors)
    return float(seen.mean()), 360.0 * longest / n_sectors


# =========================================================================== #
#  Acquisition
# =========================================================================== #

def collect_points(
    picam2,
    det,
    size,
    duration,
    centre_guess=None,
    search_half=None,
    max_step=40.0,
    fps_nominal=50.0,
):
    """
    Journalise le centre de la balle.

    La détection est effectuée dans une fenêtre carrée généreuse autour
    de centre_guess. Les coordonnées sont ensuite reconverties en
    coordonnées pleine image.
    """
    w, h = size
    xs, ys, ts = [], [], []

    if centre_guess is None or search_half is None:
        x0, y0, x1, y1 = 0, 0, w, h
        roi_desc = "image entière"
    else:
        cx0, cy0 = centre_guess
        half = float(search_half)

        x0 = max(0, int(round(cx0 - half)))
        y0 = max(0, int(round(cy0 - half)))
        x1 = min(w, int(round(cx0 + half)))
        y1 = min(h, int(round(cy0 + half)))

        roi_desc = "x=[{}, {}], y=[{}, {}]".format(x0, x1, y0, y1)

    print("[geo] zone de recherche : {}".format(roi_desc))

    t0 = time.perf_counter()
    n_frames = 0
    n_lost = 0
    n_rejected_step = 0
    next_print = 5.0

    while True:
        t = time.perf_counter() - t0
        if t >= duration:
            break

        frame = picam2.capture_array()
        n_frames += 1

        crop = frame[y0:y1, x0:x1]
        loc = det.process_image(crop)

        if loc is None:
            n_lost += 1
            continue

        # Retour aux coordonnées de l'image complète.
        px = float(loc[0]) + x0
        py = float(loc[1]) + y0

        # Filtre de vitesse / saut spatial.
        if xs:
            dt = max(t - ts[-1], 1.0 / fps_nominal)
            frame_factor = max(1.0, dt * fps_nominal)
            allowed_step = max_step * frame_factor

            jump = math.hypot(px - xs[-1], py - ys[-1])

            if jump > allowed_step:
                n_rejected_step += 1
                continue

        xs.append(px)
        ys.append(py)
        ts.append(t)

        if t >= next_print:
            print(
                "      {:4.0f} s   {:5d} points   "
                "({:d} sauts rejetés)".format(
                    t, len(xs), n_rejected_step
                )
            )
            next_print += 5.0

    print("[geo] sauts spatiaux rejetés : {}".format(n_rejected_step))

    return (
        np.asarray(xs),
        np.asarray(ys),
        np.asarray(ts),
        n_frames,
        n_lost,
        n_rejected_step,
    )

# =========================================================================== #
#  Verification appelee par les scripts d'essai
# =========================================================================== #

def check_geometry(cfg, tol_centre_px=6.0, tol_radius_px=8.0):
    """Compare la configuration au dernier fichier de calibration.

    Retourne une liste de problemes, vide si tout concorde. Ne leve rien :
    l'appelant decide si c'est bloquant, comme check_plan_against_config.

    Les tolerances sont exprimees en pixels parce que c'est l'erreur sur psi
    qui compte : 6 px sur une orbite de 208 px font 1.7 deg d'amplitude de
    modulation, du meme ordre que le bruit de detection mesure sur l'outer
    (0.2-0.6 deg) et donc deja a la limite de l'acceptable.
    """
    path = calib_path(cfg.HOOP)
    if not os.path.exists(path):
        return ["aucune calibration geometrique pour le cerceau '{}'. "
                "Lancer d'abord :\n       T1_HOOP={} python3 "
                "hoop_geometry_v3.py".format(cfg.HOOP, cfg.HOOP)]

    with open(path) as f:
        cal = json.load(f)

    problems = []
    cx_cfg, cy_cfg = cfg.HOOP_CENTRE_PX
    d = math.hypot(cx_cfg - cal["centre_px"][0], cy_cfg - cal["centre_px"][1])
    if d > tol_centre_px:
        problems.append(
            "HOOP_CENTRE_PX = ({:.1f}, {:.1f}) s'ecarte de {:.1f} px du centre "
            "mesure\n       ({:.1f}, {:.1f}) le {}. Sur une orbite de {:.0f} px "
            "cela module psi de\n       +/- {:.1f} deg le long du tour. "
            "Recopier la valeur mesuree."
            .format(cx_cfg, cy_cfg, d, cal["centre_px"][0], cal["centre_px"][1],
                    cal["timestamp"], cal["radius_px"],
                    math.degrees(d / cal["radius_px"])))

    r_cfg = float(getattr(cfg, "HOOP_RADIUS_PX", 0.0) or 0.0)
    if abs(r_cfg - cal["radius_px"]) > tol_radius_px:
        problems.append(
            "HOOP_RADIUS_PX = {:.0f} contre {:.1f} mesure : le controle de "
            "TRACKING_WIN et\n       le dimensionnement de la ROI portent sur "
            "le mauvais rayon."
            .format(r_cfg, cal["radius_px"]))

    r_roi = float(getattr(cfg, "ROI_RADIUS_PX", 0.0) or 0.0)
    besoin = cal["roi_radius_px_min"]
    if r_roi < besoin:
        problems.append(
            "ROI_RADIUS_PX = {:.0f} alors que la balle entiere en demande "
            "{:.0f} d'apres la\n       mesure geometrique : elle sera rognee "
            "au travers du cerceau.".format(r_roi, besoin))

    return problems


# =========================================================================== #
#  Programme principal
# =========================================================================== #

def main():
    from tasks import t3_config as cfg

    p = argparse.ArgumentParser(
        description="Mesure du centre et du rayon d'orbite par ajustement "
                    "de cercle")
    p.add_argument("--duration", type=float, default=30.0,
                   help="duree d'acquisition [s]")
    p.add_argument("--exposure", type=int, default=None,
                   help="pose camera [us] (defaut : celle de la config)")
    p.add_argument(
        "--search-half",
        type=float,
        default=None,
        help="demi-cote de la zone de recherche [px]",
    )

    p.add_argument(
        "--max-step",
        type=float,
        default=40.0,
        help="deplacement maximal par image [px]",
    )

    p.add_argument(
        "--k-sigma",
        type=float,
        default=3.0,
        help="seuil d'elagage robuste en ecarts-types robustes",
    )
    p.add_argument("--no-write", action="store_true",
                   help="affiche sans ecrire le fichier de calibration")
    args = p.parse_args()

    expo = args.exposure if args.exposure is not None else cfg.EXPOSURE_US

    print("=" * 70)
    print("  Geometrie du cerceau {} -- ajustement de cercle".format(cfg.HOOP))
    print("=" * 70)
    print("  image        : {} x {} px".format(*cfg.SIZE))
    print("  pose         : {} us, gain {:.1f}".format(expo, cfg.GAIN_CAM))
    print("  R_eff        : {:.5f} m (t3_config_v3.R_EFF_M)".format(cfg.R_EFF_M))
    print("  acquisition  : {:.0f} s".format(args.duration))

    centre_guess = tuple(cfg.HOOP_CENTRE_PX)

    configured_radius = float(
        getattr(cfg, "HOOP_RADIUS_PX", 270.0)
    )

    if args.search_half is None:
        search_half = max(350.0, configured_radius * 1.30)
    else:
        search_half = args.search_half

    print(
        "  centre recherche : ({:.1f}, {:.1f})".format(
            centre_guess[0], centre_guess[1]
        )
    )
    print(
        "  demi-cote recherche : {:.0f} px".format(search_half)
    )

    picam2 = open_camera(cfg.SIZE, cfg.FPS, expo, cfg.GAIN_CAM)
    # Fenetre de suivi large : le mouvement est manuel donc lent, et une
    # fenetre etroite ferait perdre la balle a la premiere hesitation.
    det = BallDetectorAA4CC(
        color_coefs=DEFAULT_COLOR_COEFS,
        threshold=cfg.THRESHOLD,
        downsample=cfg.DOWNSAMPLE,
        tracking_window=240,
        ball_size=cfg.BALL_SIZE,
    )
    for _ in range(10):
        picam2.capture_array()

    try:
        print("\n  Faites rouler la balle A LA MAIN, lentement et de facon")
        print("  reguliere, sur PLUSIEURS TOURS COMPLETS pendant toute la")
        print("  duree. Ne touchez ni la camera ni le banc.")
        input("\n  ENTREE pour demarrer ...")
        print("\n[geo] acquisition ...")
        x, y, t, n_frames, n_lost, n_step_rejected = collect_points(
            picam2,
            det,
            cfg.SIZE,
            args.duration,
            centre_guess=centre_guess,
            search_half=search_half,
            max_step=args.max_step,
            fps_nominal=cfg.FPS,
        )
    finally:
        try:
            picam2.stop()
            picam2.close()
        except Exception:
            pass

    print(
        "\n[geo] {} images, {} sans detection ({:.1f} %), "
        "{} points retenus, {} sauts rejetés".format(
            n_frames,
            n_lost,
            100.0 * n_lost / max(n_frames, 1),
            x.size,
            n_step_rejected,
        )
    )

    if x.size < 200:
        raise SystemExit(
            "trop peu de points ({}) pour un ajustement fiable. Verifier "
            "l'eclairage,\nTHRESHOLD et BALL_SIZE, puis rallonger "
            "l'acquisition.".format(x.size))

    cx, cy, R, inlier_mask, robust_history = fit_circle_robust(
        x,
        y,
        k_sigma=args.k_sigma,
        n_iter=8,
        min_keep=200,
    )

    x_fit = x[inlier_mask]
    y_fit = y[inlier_mask]

    a0, b0, R0 = fit_circle_kasa(x_fit, y_fit)
    cx, cy, R = fit_circle_geometric(
        x_fit,
        y_fit,
        (a0, b0, R0),
    )

    res = circle_residuals(x_fit, y_fit, cx, cy, R)

    rms = float(np.sqrt(np.mean(res ** 2)))
    p95 = float(np.percentile(np.abs(res), 95))
    cover, gap = angular_coverage(x_fit, y_fit, cx, cy)
    shift = math.hypot(cx - a0, cy - b0)

    n_rejected_robust = int(np.count_nonzero(~inlier_mask))

    lam = cfg.R_EFF_M / R                      # [m/px]
    r_ball_px = 0.0125 / lam
    roi_min = (R + r_ball_px) / (1.0 + cfg.ROI_MARGIN)

    print("\n  --- ajustement ---")
    print("    Kasa (algebrique)   : centre ({:7.2f}, {:7.2f})  R = {:6.2f} px"
          .format(a0, b0, R0))
    print("    Gauss-Newton (geom.): centre ({:7.2f}, {:7.2f})  R = {:6.2f} px"
          .format(cx, cy, R))
    print("    ecart entre les deux: {:.2f} px".format(shift))
    print("    residu radial       : RMS {:.2f} px, 95e centile {:.2f} px"
          .format(rms, p95))
    print("    couverture          : {:.0f} % des secteurs, plus grand trou "
          "{:.0f} deg".format(100 * cover, gap))

    print("\n --- Kasa ---")
    print(
        "    points conservés     : {} / {}".format(
            x_fit.size, x.size
        )
    )
    print(
        "    points rejetés robustes : {}".format(
            n_rejected_robust
        )
    )
    print("\n  --- grandeurs derivees ---")
    print("    lambda              : {:.4e} m/px".format(lam))
    print("    rayon de balle      : {:.1f} px (12.5 mm)".format(r_ball_px))
    print("    ROI_RADIUS_PX mini  : {:.0f} px (avec ROI_MARGIN = {:.2f})"
          .format(roi_min, cfg.ROI_MARGIN))

    # --- Diagnostics : ce qui invalide la mesure -------------------------
    fatal = []
    if gap > 60.0:
        fatal.append("un secteur de {:.0f} deg n'a jamais ete visite : le "
                     "cercle est mal contraint.".format(gap))
    if rms > 4.0:
        fatal.append("residu RMS de {:.2f} px : la balle ne decrit pas un "
                     "cercle, ou la detection derape.".format(rms))
    if shift > 5.0:
        fatal.append(
            "Kasa et Gauss-Newton different de {:.2f} px : "
            "ajustement encore instable apres rejet des aberrants."
            .format(shift)
        )
    if fatal:
        print("\n[!!] mesure NON exploitable :")
        for f in fatal:
            print("     - " + f)
        raise SystemExit("recommencer l'acquisition")

    print("\n  --- a recopier dans t1_config_v3.py / t3_config_v3.py ---")
    print("    HOOP_CENTRES_PX[{!r}] = ({:.2f}, {:.2f})".format(cfg.HOOP, cx, cy))
    print("    HOOP_RADIUS_PX        = {:.1f}".format(R))
    print("    ROI_RADIUS_PX         = {:.0f}   # >= {:.0f}"
          .format(math.ceil(roi_min / 5.0) * 5.0, roi_min))

    if args.no_write:
        return

    cal = {
        "hoop": cfg.HOOP,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "centre_px": [cx, cy],
        "radius_px": R,
        "lambda_m_per_px": lam,
        "ball_radius_px": r_ball_px,
        "roi_radius_px_min": roi_min,
        "residual_rms_px": rms,
        "residual_p95_px": p95,
        "angular_coverage": cover,
        "largest_gap_deg": gap,
        "kasa_centre_px": [a0, b0],
        "kasa_radius_px": R0,
        "n_points": int(x.size),
        "exposure_us": expo,
        "image_size": list(cfg.SIZE),
        "r_eff_m": cfg.R_EFF_M,
        "n_step_rejected": int(n_step_rejected),
        "n_robust_rejected": int(n_rejected_robust),
        "n_inliers": int(x_fit.size),
        "robust_k_sigma": float(args.k_sigma),
    }
    path = calib_path(cfg.HOOP)
    with open(path, "w") as f:
        json.dump(cal, f, indent=2)
    print("\n[out] {}".format(path))


if __name__ == "__main__":
    main()
