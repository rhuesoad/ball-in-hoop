"""
bench_common_v3.py -- Briques communes a T1 (regulation) et T2 (poursuite).

REMARQUE: traffic USB
=========================
La v2 était limitée par l'accès à l'ODrive via l'USB. Le délai était donc énorme. Attention à ce qui est demandé par cycle! 

Notations 
-----------------
    x = [theta, theta_dot, psi, psi_dot]        [rad, rad/s]
    theta = angle du cerceau, origine prise a la fermeture de boucle
    psi   = angle de la balle depuis la verticale descendante (camera)
    u     = theta_ddot                          [rad/s^2]

Loi de commande
-----------------
    u = u_ff(t) - K (x - x_ref(t))
    T1 : x_ref = [0, 0, psi_ref, 0],  u_ff = 0        (regulation)
    T2 : x_ref et u_ff issus de l'inversion du modele (poursuite)
"""

import json
import math
import os
import signal
import time
from datetime import datetime

import numpy as np
import odrive
from odrive.enums import *

from common.ball_detection import (BallDetectorAA4CC, detect_hoop_roi, open_camera,
                            DEFAULT_COLOR_COEFS)


# =========================================================================== #
#  Arret si le moindre probleme est detecte
# =========================================================================== #

class Stopper:
    """Transforme Ctrl-C en flag, pour que la boucle sorte par le chemin
    normal et beneficie de la rampe de deceleration."""

    def __init__(self):
        self.stop = False
        self.reason = ""
        signal.signal(signal.SIGINT, self._handler)
        signal.signal(signal.SIGTERM, self._handler)

    def _handler(self, signum, frame):
        self.stop = True
        self.reason = "interruption utilisateur"


# =========================================================================== #
#  Estimation de psi_dot par dérivation (arriere + passe bas ordre 1)
# =========================================================================== #

class RateEstimator:
    """
    Difference arriere et passe-bas du premier ordre (causal).
    Transmission du dt cumulé depuis la dernière détection réussie, puis dérivation.
    """

    def __init__(self, cutoff_hz):
        self.wc = 2.0 * math.pi * cutoff_hz
        self.psi_prev = None
        self.rate = 0.0

    def update(self, psi, dt):
        if self.psi_prev is None or dt <= 0.0:
            self.psi_prev = psi
            return self.rate
        raw = (psi - self.psi_prev) / dt
        self.psi_prev = psi
        a = self.wc * dt / (1.0 + self.wc * dt)   # Euler implicite : stable
        self.rate += a * (raw - self.rate)        # quel que soit dt
        return self.rate

    def reset(self):
        self.psi_prev = None
        self.rate = 0.0


# =========================================================================== #
#  Vérification de cohérence de la configuration
# =========================================================================== #

def validate_config(cfg):
    """
    Échoue s'il y a une incohérence, avant le moindre mouvement. Un réglage qui 
    se contredit lui-même doit s'arrêter avant de bouger, pas produire un essai.
    L'échec est clairement reporté. 
    """

    if cfg.COM_PERIOD != 1:
        print("\n[!] COM_PERIOD = {} : la consigne n'est renouvelee qu'a {:.1f} Hz.".format(cfg.COM_PERIOD, cfg.LOOP_HZ / cfg.COM_PERIOD))
        print("    Cela ajoute {:.0f} ms de retard equivalent.".format(500.0 * cfg.COM_PERIOD / cfg.LOOP_HZ))
        print("    Le budget USB permet COM_PERIOD = 1 : remets-le à 1 sauf si raison explicite, sinon la boucle sera dégradée :(.")

    if cfg.USE_MEASURED_THETA and cfg.ENC_PERIOD != 1:
        raise SystemExit(
            "configuration incohérente : USE_MEASURED_THETA = True exigé "
            "ENC_PERIOD = 1.\nL'encodeur entre dans la loi de commande et doit être lu a chaque tour." \
            "\n Avec ENC_PERIOD = {}, le retour d'etat voit une valeur vieille de {:.0f} ms."
            .format(cfg.ENC_PERIOD, 1e3 * cfg.ENC_PERIOD / cfg.LOOP_HZ))

    if not (0.0 < cfg.GAIN_START <= 1.0):
        raise SystemExit("GAIN_START doit etre dans ]0, 1], vaut {}."
                         .format(cfg.GAIN_START))

    if cfg.GAIN_START < 1.0 and cfg.GAIN_RAMP_S <= 0.0:
        raise SystemExit(
            "GAIN_START < 1 avec GAIN_RAMP_S = {} : le gain resterait attenue "
            "pour toujours.".format(cfg.GAIN_RAMP_S))

    if cfg.GAIN_START < 1.0:
        print("\n[!] GAIN_START = {:.2f} sur {:.2f} s : le retour d'etat est "
              "attenue".format(cfg.GAIN_START, cfg.GAIN_RAMP_S))
        print("    pendant le transitoire. C'est exactement ce qui rendait la "
              "campagne T1")
        print("    moins amortie que la simulation. A n'utiliser qu'en "
              "connaissance de cause.")

    if cfg.DT_MAX < 1.5 / cfg.LOOP_HZ:
        raise SystemExit(
            "DT_MAX = {:.3f} s est trop serre pour LOOP_HZ = {:.0f} : une "
            "seule image en retard\navorterait l'essai."
            .format(cfg.DT_MAX, cfg.LOOP_HZ))


# =========================================================================== #
#  Accès ODrive 
# =========================================================================== #

class ODriveLink:
    """Concentre TOUS les acces a l'ODrive et compte les aller-retours USB.

    Chaque lecture ou ecriture d'une propriete odrive est un aller-retour
    sur le lien. Cette communication est ce qui a le plus d'importance ici.
    Les compter ici est ce qui permet d'affirmer un budget plutot que de l'estimer.
    """

    def __init__(self, ax, cfg):
        self.ax = ax
        self.cfg = cfg
        self.n_write = 0
        self.n_read = 0
        self.t_usb = 0.0          # temps passé dans les acces ODrive [s]

        # Valeurs les plus récentes. NaN tant que rien n'a été lu : le
        # journal ne doit jamais contenir une valeur inventée.
        self.theta_enc = float("nan")
        self.theta_dot_enc = float("nan")
        self.iq = float("nan")
        self.errors = (0, 0)

    # -- Ecriture ---------------------
    def send_velocity(self, theta_dot_cmd):
        t0 = time.perf_counter()
        self.ax.controller.input_vel = theta_dot_cmd / (2 * math.pi)
        self.t_usb += time.perf_counter() - t0
        self.n_write += 1

    # -- Lectures ------------
    def read_encoder(self, theta0):
        """Retourne (theta, theta_dot) en rad et rad/s, ou (nan, nan) si ce
        tour n'est pas un tour de lecture. L'appelant decide."""
        t0 = time.perf_counter()
        pos = self.ax.pos_estimate
        vel = self.ax.vel_estimate
        self.t_usb += time.perf_counter() - t0
        self.n_read += 2
        self.theta_enc = (pos - theta0) * 2 * math.pi
        self.theta_dot_enc = vel * 2 * math.pi
        return self.theta_enc, self.theta_dot_enc

    def read_iq(self):
        t0 = time.perf_counter()
        self.iq = self.ax.motor.foc.Iq_measured
        self.t_usb += time.perf_counter() - t0
        self.n_read += 1
        return self.iq

    def read_errors(self):
        t0 = time.perf_counter()
        active = self.ax.active_errors
        disarm = self.ax.disarm_reason
        self.t_usb += time.perf_counter() - t0
        self.n_read += 2
        self.errors = (active, disarm)
        return self.errors

    def read_position(self):
        """Lecture ponctuelle de la position, pour fixer l'origine theta0."""
        t0 = time.perf_counter()
        pos = self.ax.pos_estimate
        self.t_usb += time.perf_counter() - t0
        self.n_read += 1
        return pos

    # -- Budget ------------------------------------------------------------
    def report(self, n_ticks):
        if n_ticks <= 0:
            return {}
        n_ops = self.n_write + self.n_read
        return {
            "usb_writes": self.n_write,
            "usb_reads": self.n_read,
            "usb_ops_per_tick": n_ops / n_ticks,
            "usb_ms_per_tick": 1e3 * self.t_usb / n_ticks,
        }


# =========================================================================== #
#  ODrive : connexion et arret
# =========================================================================== #

def preflight_check(odrv0, ax, cfg):
    """Verifications en lecture seule avant tout mouvement. Hors boucle,
    donc le trafic USB n'y est pas un enjeu."""
    br = odrv0.config.brake_resistor0
    if not br.enable:
        raise RuntimeError(
            "resistance de freinage desactivee : le bus DC montera en tension "
            "au freinage. Activer via odrivetool puis save_configuration().")
    print("[odrv] frein : {:.1f} ohm, active".format(br.resistance))

    for path in ("fet_thermistor", "motor_thermistor"):
        th = getattr(ax.motor, path, None)
        if th is None:
            continue
        cfg_th = getattr(th, "config", None)
        print("[odrv] {} : {:.1f} degC (enabled={})".format(
            path, getattr(th, "temperature", float("nan")),
            getattr(cfg_th, "enabled", "?")))

    print("[odrv] Vbus = {:.2f} V".format(odrv0.vbus_voltage))
    print("[odrv] limite de courant : {:.1f} A  (soit {:.3f} N.m a Kt={:.3f})"
          .format(ax.config.motor.current_soft_max,
                  ax.config.motor.current_soft_max * cfg.KT, cfg.KT))


def connect_and_prepare(cfg):
    print("[odrv] connexion ...")
    odrv0 = odrive.find_any()
    ax = odrv0.axis0
    print("[odrv] connecté (serial {}, fw {}.{}.{})".format(
        odrv0.serial_number, odrv0.fw_version_major,
        odrv0.fw_version_minor, odrv0.fw_version_revision))

    preflight_check(odrv0, ax, cfg)
    # L'AMT102-V n'a pas d'index : l'offset encodeur est perdu a chaque
    # coupure du driver. La calibration moteur (R, L), elle, est persistante.
    if ax.current_state != AXIS_STATE_CLOSED_LOOP_CONTROL:
        print("[odrv] calibration de l'offset encodeur ...")
        ax.requested_state = AXIS_STATE_ENCODER_OFFSET_CALIBRATION
        while ax.current_state != AXIS_STATE_IDLE:
            time.sleep(0.1)
        if ax.active_errors != 0 or ax.disarm_reason != 0:
            raise RuntimeError(
                "calibration echouee : active=0x{:X}, disarm=0x{:X}"
                .format(ax.active_errors, ax.disarm_reason))
        ax.requested_state = AXIS_STATE_CLOSED_LOOP_CONTROL
        time.sleep(0.5)

    ax.controller.config.control_mode = CONTROL_MODE_VELOCITY_CONTROL
    ax.controller.config.input_mode = INPUT_MODE_PASSTHROUGH

    # Valeur IMPOSÉE, sinon la limite depend de ce qui traine dans la configuration
    # auvegardée et change d'un essai a l'autre. C'est elle qui fixe l'enveloppe de 
    # poursuite de T2.

    vel_limit_turns = cfg.THETA_DOT_MAX / (2 * math.pi)
    ax.controller.config.vel_limit = 1.2 * vel_limit_turns
    print("[odrv] vel_limit impose a {:.2f} tr/s (soit {:.1f} rad/s), "
          "consigne plafonnee a {:.1f} rad/s"
          .format(1.2 * vel_limit_turns, 1.2 * cfg.THETA_DOT_MAX,
                  cfg.THETA_DOT_MAX))

    ax.controller.input_vel = 0.0

    if ax.active_errors != 0 or ax.disarm_reason != 0:
        raise RuntimeError("ODrive en erreur avant demarrage : "
                           "active=0x{:X}, disarm=0x{:X}"
                           .format(ax.active_errors, ax.disarm_reason))

    print("[odrv] mode vitesse arme, consigne a zero")
    return odrv0, ax


def ramp_down(ax, theta_dot_cmd, cfg):
    """Deceleration progressive puis mise en IDLE.

    Couper la commande d'un coup laisserait le cerceau lance : l'ODrive
    encaisserait le freinage en butee de courant, ce qui peut declencher
    CURRENT_LIMIT_VIOLATION et laisser l'axe désarmé dans un état sale.
    """
    n = max(int(cfg.RAMP_DOWN_S * cfg.LOOP_HZ), 1)
    for i in range(n):
        v = theta_dot_cmd * (1.0 - (i + 1) / n)
        try:
            ax.controller.input_vel = v / (2 * math.pi)
        except Exception:
            break
        time.sleep(cfg.DT_NOMINAL)
    try:
        ax.controller.input_vel = 0.0
        time.sleep(0.2)
        if not getattr(cfg, "KEEP_ARMED", False):
            ax.requested_state = AXIS_STATE_IDLE
    except Exception:
        pass


# =========================================================================== #
#  Camera
# =========================================================================== #

def roi_from_configured_hoop(cfg):
    """Construit une ROI autour du hoop choisi.
    """
    cx, cy = cfg.HOOP_CENTRE_PX
    radius = getattr(cfg, "ROI_RADIUS_PX", cfg.HOOP_RADIUS_PX)
    margin = cfg.ROI_MARGIN

    image_w, image_h = cfg.SIZE
    half = radius * (1.0 + margin)

    x0 = max(0, int(round(cx - half)))
    y0 = max(0, int(round(cy - half)))
    x1 = min(image_w, int(round(cx + half)))
    y1 = min(image_h, int(round(cy + half)))

    return (
        x0,
        y0,
        x1 - x0,
        y1 - y0,
    )

def setup_camera(cfg):
    picam2 = open_camera(cfg.SIZE, cfg.FPS, cfg.EXPOSURE_US, cfg.GAIN_CAM)
    det = BallDetectorAA4CC(
        color_coefs=DEFAULT_COLOR_COEFS,
        threshold=cfg.THRESHOLD,
        downsample=cfg.DOWNSAMPLE,
        tracking_window=cfg.TRACKING_WIN,
        ball_size=cfg.BALL_SIZE,
    )
    for _ in range(10):
        picam2.capture_array()

    # Le ROI issu de Hough ne sert qu'a recadrer la recherche : sa precision
    # est sans consequence. Seul le centre entre dans le calcul de psi.
    roi = detect_hoop_roi(picam2.capture_array(), margin=cfg.ROI_MARGIN)
    rx, ry, rw, rh = roi

    if cfg.USE_HOUGH_CENTRE:
        roi = roi_from_configured_hoop(cfg)
        centre = tuple(float(v) for v in cfg.HOOP_CENTRE_PX)

        print(
            "[cam] hoop={} centre configure=({:.1f}, {:.1f}) px".format(
                cfg.HOOP,
                centre[0],
                centre[1],
            )
        )
    else:
        roi = detect_hoop_roi(
            picam2.capture_array(),
            margin=cfg.ROI_MARGIN,
        )
        rx, ry, rw, rh = roi
        centre = (
            rx + rw / 2.0,
            ry + rh / 2.0,
        )

    rx, ry, rw, rh = roi

    print(
        "[cam] ROI hoop={} : x[{}:{}] y[{}:{}] ({} x {} px)".format(
            cfg.HOOP,
            rx,
            rx + rw,
            ry,
            ry + rh,
            rw,
            rh,
        )
    )

    print(
        "[cam] centre utilisé pour psi : "
        "({:.1f}, {:.1f}) px".format(
            centre[0],
            centre[1],
        )
    )

    return picam2, det, roi, centre


def _read_psi_raw(picam2, det, roi, centre):
    """Une image -> psi BRUT (offset non retranche), ou None."""
    rx, ry, rw, rh = roi
    cx, cy = centre
    frame = picam2.capture_array()
    loc = det.process_image(frame[ry:ry + rh, rx:rx + rw])
    if loc is None:
        return None
    u_px, v_px = loc[0] + rx, loc[1] + ry
    return math.atan2(u_px - cx, v_px - cy)


# =========================================================================== #
#  Armement
# =========================================================================== #

class ArmingRetry(Exception):
    """Condition d'armement non satisfaite, mais rattrapable par
    l'operateur : on repropose au lieu d'abandonner l'essai."""


def _retry_prompt(exc, attempt=[0]):
    attempt[0] += 1
    print("\n[arm] CONDITION NON REMPLIE (tentative {}) :".format(attempt[0]))
    for line in str(exc).splitlines():
        print("      " + line)
    ans = input("\n      ENTREE pour recommencer, ou 'q' + ENTREE pour "
                "abandonner : ").strip().lower()
    return ans != "q"


def psi_statistics(t, psi, f_n_hz):
    """
    Sépare mouvement et bruit dans une série de psi mesurés au repos.

    Une balle qui bouge encore oscille à la fréquence propre f_n. 
    Le bruit de détection est large bande. 
    Une régression linéaire permet de les séparer. 

    psi(t) ~= a sin(w_n t) + b cos(w_n t) + c

    - amplitude sqrt(a^2 + b^2)  -> mouvement residuel reel
    - c                          -> offset camera, non biaise par
                                    l'oscillation puisqu'elle est
                                    explicitement modelisee
    - ecart-type du résidu       -> bruit de detection, insensible
                                    aux detections aberrantes

    Dériver psi serait mauvais, car diviser par dt = 20 ms amplifie 
    le bruit de mesure d'un facteur 50, donc RMS(dpsi/dt) vaudrait 
    1.24 rad/s pour sigma_psi = 1 deg, balle parfaitement immobile.

    [Ljung, "System Identification: Theory for the User", 2e ed., 1999]
    """
    
    t = np.asarray(t, float)
    psi = np.asarray(psi, float)
    w = 2.0 * math.pi * f_n_hz

    M = np.column_stack([np.sin(w * t), np.cos(w * t), np.ones(t.size)])
    coef, *_ = np.linalg.lstsq(M, psi, rcond=None)
    resid = psi - M @ coef

    mad = float(np.median(np.abs(resid - np.median(resid))))
    sigma = 1.4826 * mad
    outliers = float(np.mean(np.abs(resid) > 5.0 * sigma)) if sigma > 0 else 0.0

    return {
        "offset": float(coef[2]),
        "median": float(np.median(psi)),
        "osc_amp": float(math.hypot(coef[0], coef[1])),
        "sigma": float(sigma),
        "p2p": float(np.max(psi) - np.min(psi)),
        "outlier_frac": outliers,
        "n": int(t.size),
    }


def _observe(picam2, det, roi, centre, duration):
    """Collecte psi brut pendant duration. Retourne (t, psi) en tableaux."""
    samples, times = [], []
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < duration:
        psi_raw = _read_psi_raw(picam2, det, roi, centre)
        t = time.perf_counter()
        if psi_raw is None:
            continue
        samples.append(psi_raw)
        times.append(t - t0)

    if len(samples) < 20:
        raise ArmingRetry(
            "balle detectee sur {} images seulement.\n"
            "Verifier qu'elle est bien dans le champ, l'eclairage et le\n"
            "seuil de detection (THRESHOLD).".format(len(samples)))
    return np.asarray(times), np.asarray(samples)


def _arm_step1(picam2, det, roi, centre, cfg):
    """Mesure de l'offset camera, balle au repos au fond."""
    print("\n[arm] ETAPE 1 : posez la balle AU FOND du cerceau et laissez-la "
          "s'immobiliser.")
    input("      ENTREE quand elle est au repos ...")

    print("[arm] observation pendant {:.1f} s ...".format(cfg.ARM_SETTLE_S))
    t, psi = _observe(picam2, det, roi, centre, cfg.ARM_SETTLE_S)
    st = psi_statistics(t, psi, cfg.F_N_HZ)

    print("[arm] {} images".format(st["n"]))
    print("      offset camera        : {:+6.2f} deg   (limite {:.1f})"
          .format(math.degrees(st["offset"]), cfg.ARM_PSI_MAX_DEG))
    print("      oscillation a {:.2f} Hz : {:6.2f} deg   (limite {:.1f})"
          .format(cfg.F_N_HZ, math.degrees(st["osc_amp"]),
                  cfg.ARM_OSC_MAX_DEG))
    print("      bruit de detection   : {:6.2f} deg   (limite {:.1f})"
          .format(math.degrees(st["sigma"]), cfg.ARM_NOISE_MAX_DEG))
    print("      crete-a-crete        : {:6.2f} deg,  aberrants {:.1f} %"
          .format(math.degrees(st["p2p"]), 100 * st["outlier_frac"]))

    if math.degrees(st["osc_amp"]) > cfg.ARM_OSC_MAX_DEG:
        raise ArmingRetry(
            "la balle oscille encore : {:.2f} deg d'amplitude a {:.2f} Hz\n"
            "(limite {:.1f} deg). Attendre l'immobilisation."
            .format(math.degrees(st["osc_amp"]), cfg.F_N_HZ,
                    cfg.ARM_OSC_MAX_DEG))

    if math.degrees(st["sigma"]) > cfg.ARM_NOISE_MAX_DEG:
        raise ArmingRetry(
            "detection trop bruitee : {:.2f} deg (limite {:.1f}).\n"
            "Ce n'est PAS la balle : c'est la vision. Verifier l'eclairage,\n"
            "les reflets sur le cerceau, THRESHOLD et BALL_SIZE. Lancer\n"
            "psi_bench_v3.py pour caracteriser."
            .format(math.degrees(st["sigma"]), cfg.ARM_NOISE_MAX_DEG))

    if abs(math.degrees(st["offset"])) > cfg.ARM_PSI_MAX_DEG:
        raise ArmingRetry(
            "offset camera de {:+.1f} deg, au-dela de {:.1f} deg.\n"
            "La physique impose psi = 0 pour une balle au repos au fond :\n"
            "soit la balle n'est pas au fond, soit le centre du cerceau en\n"
            "pixels (HOOP_CENTRE_PX) est faux."
            .format(math.degrees(st["offset"]), cfg.ARM_PSI_MAX_DEG))

    return st["offset"], st


def _arm_step2(picam2, det, roi, centre, cfg, psi_offset, psi_init_deg):
    """Mise en place de la condition initiale.
    """
    print("\n[arm] ETAPE 2 : ecartez la balle a {:+.0f} deg et MAINTENEZ-LA."
          .format(psi_init_deg))
    input("      ENTREE quand elle est en place ...")

    t, psi = _observe(picam2, det, roi, centre, cfg.ARM_HOLD_S)
    psi_now = float(np.median(psi)) - psi_offset
    psi_now_deg = math.degrees(psi_now)
    print("[arm] psi mesure : {:+.1f} deg (vise {:+.1f}, ecart {:+.1f})"
          .format(psi_now_deg, psi_init_deg, psi_now_deg - psi_init_deg))

    if abs(psi_now_deg - psi_init_deg) > cfg.ARM_INIT_TOL_DEG:
        raise ArmingRetry(
            "ecart de {:+.1f} deg a la consigne initiale (tolerance "
            "{:.1f} deg).\nRepositionner la balle."
            .format(psi_now_deg - psi_init_deg, cfg.ARM_INIT_TOL_DEG))
    return psi_now_deg


def arm(picam2, det, roi, centre, cfg, psi_init_deg):
    """Mesure l'offset camera puis amene le banc a la condition initiale.

    Retourne (psi_offset, psi_init_mesure_deg, diagnostics_etape1).
    """
    while True:
        try:
            psi_offset, stats = _arm_step1(picam2, det, roi, centre, cfg)
            break
        except ArmingRetry as e:
            if not _retry_prompt(e):
                raise RuntimeError("armement abandonne par l'operateur.")

    if abs(psi_init_deg) < 1e-6:
        print("\n[arm] ETAPE 2 : condition initiale psi = 0, rien a faire.")
        input("      ENTREE pour fermer la boucle (Ctrl-C pour arreter) ...")
        return psi_offset, 0.0, stats

    while True:
        try:
            psi_init_meas = _arm_step2(picam2, det, roi, centre, cfg,
                                       psi_offset, psi_init_deg)
            break
        except ArmingRetry as e:
            if not _retry_prompt(e):
                raise RuntimeError("armement abandonne par l'operateur.")

    for c in (3, 2, 1):
        print("      {} ...".format(c))
        time.sleep(1.0)
    print("      LACHEZ")
    return psi_offset, psi_init_meas, stats


# =========================================================================== #
#  Boucle de commande
# =========================================================================== #

def control_loop(ax, picam2, det, roi, centre, psi_offset,
                 K, reference, duration, stopper, cfg,
                 estimator=None, tag="t1", gain_schedule=None,
                 phase_hook=None):
    """Boucle a LOOP_HZ, cadencee par la camera.
    """
    rx, ry, rw, rh = roi
    cx, cy = centre

    rate_est = RateEstimator(cfg.PSIDOT_FILTER_HZ)
    link = ODriveLink(ax, cfg)

    theta0 = link.read_position()            # origine de theta, hors boucle

    unwrap_psi = bool(getattr(cfg, "UNWRAP_PSI", False))
    psi_acc = None                  # None tant qu'aucune image n'est vue
    psi_raw_prev = 0.0

    theta_int = 0.0                 # integrale de theta_dot_hat
    theta_dot_cmd = 0.0             # integrale de u  (= consigne de vitesse)

   
    tau_vel = float(getattr(cfg, "TAU_VEL_S", 0.0) or 0.0)
    theta_dot_hat = 0.0

    lost = 0
    overrun = 0
    dt_meas = 0.0                   # dt cumule depuis la derniere mesure
    stop_reason = "duree atteinte"

    n_max = int(duration * cfg.LOOP_HZ * 1.5) + 200
    keys = ("t", "dt", "theta", "theta_dot", "theta_enc", "theta_dot_enc",
            "psi", "psi_dot", "psi_ref", "psi_dot_ref", "theta_ref",
            "theta_dot_ref", "u_ff", "u_raw", "u", "vel_cmd", "iq",
            "gain_scale", "n_lost", "r_px", "phase")
    log = {k: np.full(n_max, np.nan) for k in keys}
    i = 0

    seg = {k: 0.0 for k in ("wait", "det", "ctrl", "io")}
    seg_n = 0

    COM_PERIOD = cfg.COM_PERIOD
    ENC_PERIOD = cfg.ENC_PERIOD
    IQ_PERIOD = cfg.IQ_PERIOD
    ERR_PERIOD = cfg.ERR_PERIOD

    print("\n[{}] boucle fermee -- Ctrl-C pour arreter\n".format(tag))

    t_start = time.perf_counter()
    t_prev = t_start

    try:
        while True:
            # ---- Camera (bloquante : c'est elle qui cadence la boucle) ----
            t_a = time.perf_counter()
            frame = picam2.capture_array()
            t_b = time.perf_counter()

            dt = t_b - t_prev
            t_prev = t_b
            t_rel = t_b - t_start
            dt_meas += dt

            if stopper.stop:
                stop_reason = stopper.reason
                break
            if t_rel >= duration or i >= n_max:
                break

            if dt > cfg.DT_MAX:
                overrun += 1
                if overrun >= cfg.ABORT_OVERRUN:
                    stop_reason = ("budget de boucle depasse {} fois de suite "
                                   "(dernier dt = {:.0f} ms)"
                                   .format(overrun, dt * 1e3))
                    break
            else:
                overrun = 0

            # ---- Detection ------------------------------------------------
            loc = det.process_image(frame[ry:ry + rh, rx:rx + rw])
            t_c = time.perf_counter()

            if loc is None:
                lost += 1
                if lost >= cfg.ABORT_LOST_FRAMES:
                    stop_reason = ("balle perdue pendant {} images "
                                   "consecutives".format(lost))
                    break

                if estimator is not None:
                    estimator.predict_open_loop(dt_meas)
                    estimator.push_command(0.0)
                    dt_meas = 0.0
                if (i % COM_PERIOD) == 0:
                    link.send_velocity(theta_dot_cmd)
                continue
            n_lost_here = lost
            lost = 0

            u_px, v_px = loc[0] + rx, loc[1] + ry
            psi_raw = math.atan2(u_px - cx, v_px - cy) - psi_offset
            r_px = math.hypot(u_px - cx, v_px - cy)

            if unwrap_psi:
                if psi_acc is None:
                    psi_acc = psi_raw
                else:
                    d = psi_raw - psi_raw_prev
                    d -= 2.0 * math.pi * round(d / (2.0 * math.pi))
                    psi_acc += d
                psi_raw_prev = psi_raw
                psi = psi_acc
            else:
                psi = psi_raw

            if abs(psi) > math.radians(cfg.ABORT_PSI_DEG):
                stop_reason = ("psi = {:+.1f} deg, au-dela de la limite de "
                               "{:.0f} deg".format(math.degrees(psi),
                                                   cfg.ABORT_PSI_DEG))
                break

            # ---- Estimation de psi_dot ------------------------------------
            if estimator is None:
                psi_dot = rate_est.update(psi, dt_meas)
            else:
                xk = estimator.update(psi, dt_meas)
                psi, psi_dot = float(xk[0]), float(xk[1])
            dt_meas = 0.0

            # ---- Loi de commande ------------------------------------------

            if cfg.GAIN_START >= 1.0:
                gain_scale = 1.0
            else:
                gain_scale = min(1.0, cfg.GAIN_START +
                                 (1.0 - cfg.GAIN_START) * t_rel / cfg.GAIN_RAMP_S)

            if cfg.USE_MEASURED_THETA:
                theta, theta_dot = link.theta_enc, link.theta_dot_enc
            else:
                theta, theta_dot = theta_int, theta_dot_hat

            x = np.array([theta, theta_dot, psi, psi_dot])

            phase = 0
            if phase_hook is None:
                x_ref, u_ff = reference(t_rel)
                K_t = K if gain_schedule is None else gain_schedule(t_rel)
            else:
                x_ref, u_ff, K_t, phase = phase_hook(t_rel, x, r_px)

            if K_t is None:
                u_raw = u_ff              # boucle ouverte assumee (T4 phase 2)
            else:
                u_raw = u_ff - gain_scale * float(K_t @ (x - x_ref))

            # ---- Saturations ----------------------------------------------
            # 1) couple disponible (asymetrique : le frottement visqueux aide
            #    dans un sens et s'oppose dans l'autre)
            u_hi = (cfg.TAU_MAX - cfg.B_MOTEUR * theta_dot) / cfg.I_TOTAL
            u_lo = (-cfg.TAU_MAX - cfg.B_MOTEUR * theta_dot) / cfg.I_TOTAL
            u = min(max(u_raw, u_lo), u_hi)
            # 2) borne de non-glissement balle / joints toriques
            u = min(max(u, -cfg.U_SLIP_MAX), cfg.U_SLIP_MAX)

            # ---- Integration en consigne de vitesse (anti-windup) ---------
            v_new = theta_dot_cmd + u * dt
            if abs(v_new) > cfg.THETA_DOT_MAX:
                v_new = math.copysign(cfg.THETA_DOT_MAX, v_new)
                u = (v_new - theta_dot_cmd) / dt     # back-calculation
            # La consigne envoyee au variateur reste v_new : c'est elle qui
            # ferme la boucle. Seule l'ESTIMATION servant au retour d'etat
            # est retardee du premier ordre identifie.
            theta_dot_cmd = v_new
            if tau_vel > 0.0:
                # Euler implicite : stable pour tout dt, comme RateEstimator.
                a_vel = dt / (tau_vel + dt)
                hat_new = theta_dot_hat + a_vel * (theta_dot_cmd - theta_dot_hat)
            else:
                hat_new = theta_dot_cmd
            # trapeze : coherent avec theta_ddot constant sur le pas
            theta_int += 0.5 * (theta_dot_hat + hat_new) * dt
            theta_dot_hat = hat_new

            t_d = time.perf_counter()

            # ---- Acces ODrive ---------------------------------------------
            # La consigne d'abord, et a chaque tour : c'est le seul acces qui
            # ferme la boucle. Tout le reste est de la surveillance.
            if (i % COM_PERIOD) == 0:
                link.send_velocity(theta_dot_cmd)

            enc_fresh = (i % ENC_PERIOD) == 0
            if enc_fresh:
                link.read_encoder(theta0)

            if IQ_PERIOD > 0 and (i % IQ_PERIOD) == 0:
                link.read_iq()

            if (i % ERR_PERIOD) == 0:
                active, disarm = link.read_errors()
                if active != 0 or disarm != 0:
                    stop_reason = ("erreur ODrive : active=0x{:X}, "
                                   "disarm=0x{:X}".format(active, disarm))
                    break

            if estimator is not None:
                estimator.push_command(u)
            for key, val in (
                    ("t", t_rel), ("dt", dt),
                    ("theta", theta), ("theta_dot", theta_dot),
                    ("psi", psi), ("psi_dot", psi_dot),
                    ("theta_ref", x_ref[0]), ("theta_dot_ref", x_ref[1]),
                    ("psi_ref", x_ref[2]), ("psi_dot_ref", x_ref[3]),
                    ("u_ff", u_ff), ("u_raw", u_raw), ("u", u),
                    ("vel_cmd", theta_dot_cmd),
                    ("gain_scale", gain_scale), ("n_lost", n_lost_here),
                    ("r_px", r_px), ("phase", phase)):
                log[key][i] = val
            if enc_fresh:
                log["theta_enc"][i] = link.theta_enc
                log["theta_dot_enc"][i] = link.theta_dot_enc
            if IQ_PERIOD > 0 and (i % IQ_PERIOD) == 0:
                log["iq"][i] = link.iq
            i += 1

            t_e = time.perf_counter()

            seg["wait"] += t_b - t_a
            seg["det"] += t_c - t_b
            seg["ctrl"] += t_d - t_c
            seg["io"] += t_e - t_d
            seg_n += 1

            if i % int(cfg.LOOP_HZ) == 0:
                print("\r  t={:5.1f}s  psi={:+6.1f}  ref={:+6.1f}  "
                      "err={:+6.2f} deg  u={:+7.1f}  vel={:+6.2f} rad/s "
                      .format(t_rel, math.degrees(psi),
                              math.degrees(x_ref[2]),
                              math.degrees(psi - x_ref[2]),
                              u, theta_dot_cmd), end="", flush=True)

    finally:
        print("\n[{}] arret : {}".format(tag, stop_reason))
        ramp_down(ax, theta_dot_cmd, cfg)
        print("[{}] cerceau arrete, ODrive en IDLE".format(tag))

        if seg_n > 0:
            tot = sum(seg.values()) / seg_n * 1e3
            print("\n  budget de boucle sur {} pas (segments disjoints) :"
                  .format(seg_n))
            for k, label in (("wait", "attente image (capture)"),
                             ("det", "detection balle"),
                             ("ctrl", "loi de commande"),
                             ("io", "acces ODrive + journal")):
                ms = seg[k] / seg_n * 1e3
                print("    {:26s} {:6.2f} ms   ({:4.1f} %)"
                      .format(label, ms, 100.0 * ms / tot if tot else 0.0))
            print("    {:26s} {:6.2f} ms   -> {:.1f} Hz"
                  .format("TOTAL", tot, 1e3 / tot if tot else 0.0))

            usb = link.report(seg_n)
            print("\n  budget USB :")
            print("    ecritures (consigne)       {:6d}".format(usb["usb_writes"]))
            print("    lectures  (enc/Iq/erreurs) {:6d}".format(usb["usb_reads"]))
            print("    aller-retours par tour     {:6.2f}   "
                  "(v2 : 3.70)".format(usb["usb_ops_per_tick"]))
            print("    temps USB par tour         {:6.2f} ms  "
                  "({:.1f} % de la periode)"
                  .format(usb["usb_ms_per_tick"],
                          100.0 * usb["usb_ms_per_tick"] * cfg.LOOP_HZ / 1e3))
            log["_usb"] = usb

    for key in list(log):
        if key.startswith("_"):
            continue
        log[key] = log[key][:i]
    log["stop_reason"] = stop_reason
    return log


# =========================================================================== #
#  Sauvegarde et resume
# =========================================================================== #

def save_and_summarise(log, cfg, outdir, name_fmt, meta_extra, K, gain_name):
    os.makedirs(outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(outdir, name_fmt.format(stamp=stamp))

    usb = log.pop("_usb", {})

    meta = {
        "timestamp": stamp,
        "bench_version": 3,
        "gain_name": gain_name,
        "K": list(map(float, K)),
        "use_measured_theta": bool(cfg.USE_MEASURED_THETA),
        "gain_start": cfg.GAIN_START, "gain_ramp_s": cfg.GAIN_RAMP_S,
        "tau_max": cfg.TAU_MAX, "b_moteur": cfg.B_MOTEUR,
        "I_total": cfg.I_TOTAL, "kt": cfg.KT,
        "u_slip_max": cfg.U_SLIP_MAX, "theta_dot_max": cfg.THETA_DOT_MAX,
        "psidot_filter_hz": cfg.PSIDOT_FILTER_HZ,
        "cam_latency_s": cfg.CAM_LATENCY_S,
        "loop_hz": cfg.LOOP_HZ,
        "com_period": cfg.COM_PERIOD, "enc_period": cfg.ENC_PERIOD,
        "iq_period": cfg.IQ_PERIOD, "err_period": cfg.ERR_PERIOD,
        "a21": cfg.A21, "a22": cfg.A22, "b2": cfg.B2,
        "hoop_centre_px": list(cfg.HOOP_CENTRE_PX),
        "stop_reason": log.pop("stop_reason"),
        "tau_vel_s": float(getattr(cfg, "TAU_VEL_S", 0.0) or 0.0),
    }
    meta.update(usb)
    meta.update(meta_extra)
    np.savez_compressed(path, meta=json.dumps(meta), **log)

    n = log["t"].size
    if n > int(cfg.LOOP_HZ):
        tail = slice(int(n * 0.5), n)
        err = np.degrees(log["psi"][tail] - log["psi_ref"][tail])
        rate = 1.0 / np.median(log["dt"][1:n])
        print("\n  echantillons        : {}".format(n))
        print("  cadence mediane     : {:.1f} Hz   (visee {:.0f})"
              .format(rate, cfg.LOOP_HZ))
        if rate < 0.95 * cfg.LOOP_HZ:
            print("  [!] la boucle n'a pas tenu la cadence. Regarder le budget")
            print("      ci-dessus : 'attente image' dominant = camera, ")
            print("      'acces ODrive' dominant = lien USB.")
        # Pour une reference a moyenne nulle, la moyenne de l'erreur ne dit
        # rien : c'est le RMS qui compte.
        print("  erreur (2e moitie)  : RMS {:.2f} deg, crete {:.2f} deg, "
              "biais {:+.2f} deg"
              .format(float(np.sqrt(np.mean(err ** 2))),
                      float(np.max(np.abs(err))), float(err.mean())))
        print("  |psi| max           : {:.1f} deg"
              .format(np.degrees(np.abs(log["psi"])).max()))
        print("  |u| max             : {:.1f} rad/s^2  (borne {:.0f})"
              .format(np.abs(log["u"]).max(), cfg.U_SLIP_MAX))
        print("  saturation de u     : {:.1f} % des pas"
              .format(100.0 * np.mean(np.abs(log["u_raw"] - log["u"]) > 1e-6)))
        print("  |vel_cmd| max       : {:.2f} rad/s  (borne {:.0f})"
              .format(np.abs(log["vel_cmd"]).max(), cfg.THETA_DOT_MAX))
        # Ecart entre la consigne de vitesse et la vitesse reelle : mesure
        # directe de la validite de l'hypothese "boucle interne parfaite".
        # Sur les tours sans lecture encodeur les deux termes sont NaN, d'ou
        # les reductions nan-safe.
        d = log["vel_cmd"] - log["theta_dot_enc"]
        if np.any(np.isfinite(d)):
            print("  |vel_cmd - vel_enc| : RMS {:.3f} rad/s, max {:.3f} rad/s "
                  "(sur {} lectures)"
                  .format(float(np.sqrt(np.nanmean(d ** 2))),
                          float(np.nanmax(np.abs(d))),
                          int(np.sum(np.isfinite(d)))))
        if np.any(np.isfinite(log["iq"])):
            print("  |Iq| max            : {:.2f} A"
                  .format(float(np.nanmax(np.abs(log["iq"])))))
        print("  images perdues      : {:.0f}".format(np.nansum(log["n_lost"])))
    print("\n[out] {}".format(path))
    return path


def shutdown(ax, picam2, cfg=None):
    """cfg optionnel : sans lui, comportement inchange (desarmement)."""
    if ax is not None:
        try:
            ax.controller.input_vel = 0.0
            if cfg is None or not getattr(cfg, "KEEP_ARMED", False):
                ax.requested_state = AXIS_STATE_IDLE
        except Exception:
            pass
    try:
        picam2.stop()
        picam2.close()
    except Exception:
        pass


def pick_gain(cfg, name):
    """Un seul point d'entree pour le choix des gains : le nom passe en
    ligne de commande est la cle du dictionnaire, il est affiche et il est
    enregistre dans meta."""
    if name not in cfg.GAINS:
        raise SystemExit("jeu de gains inconnu : {} (disponibles : {})"
                         .format(name, ", ".join(sorted(cfg.GAINS))))
    return np.asarray(cfg.GAINS[name], dtype=float)
