#!/usr/bin/env python3
"""
e4_run.py -- E4 : acquisition de la campagne de mesure de mu.

    python3 e4_run.py --dry-run           # verifie la config, ne bouge pas
    python3 e4_run.py                     # 10 essais : 5 en +, 5 en -
    python3 e4_run.py --direction +1 --trials 2
    python3 e4_run.py --resume 7           # reprend a l'essai 7

DEROULEMENT D'UN ESSAI
======================
1. Retour a theta = 0, bille laissee redescendre au fond.
2. Mesure de l'offset camera, bille au repos : la physique impose psi = 0,
   donc ce qui est mesure est l'erreur d'alignement du centre. Meme
   procedure et memes criteres d'acceptation que bc.arm() etape 1.
3. Boucle d'increment : theta += direction * STEP_DEG, attente de
   stabilisation de la bille, enregistrement de (theta, psi) sur RECORD_S.
4. A partir du pas N_FIT_STEPS, ajustement de psi = a*theta + b sur les
   premiers pas et surveillance du residu. Deux pas consecutifs au-dela de
   SLIP_RESID_DEG avec residu croissant = decrochage ; theta* est l'angle
   du dernier pas d'avant.
5. Quelques pas de plus (pour le graphe), puis retour a zero.

CE QUI EST JOURNALISE
=====================
Le flux CONTINU (t, theta, psi, indice de pas, phase) a 50 Hz pendant tout
l'essai, ET les agregats par pas. Le flux continu coute quelques milliers
de points et il est le seul moyen de rejuger a posteriori un critere de
stabilisation ou de decrochage sans refaire la manip.

CE QUE CE SCRIPT NE FAIT PAS
============================
Il ne reconfigure NI la limite de courant NI la calibration de l'ODrive.
Il verifie que la limite en place est celle des campagnes precedentes et
refuse de partir sinon. Un script d'essai qui ecrit dans motor.config
laisse le banc dans un etat different de celui ou les autres essais ont
ete faits, et cela ne se voit dans aucun journal.
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

from experiments import e4_config as cfg
from common import bench_common as bc


# ============================================================================
# VERIFICATIONS PREALABLES
# ============================================================================

def check_geometry_consistency():
    """L'annexe A et les cotes as-built doivent donner le meme contact.

    MU_K et MU_ALPHA_DEG sont les valeurs CITEES dans le memoire ; elles
    sont recalculees dans e4_config a partir de Rb, rc et h. Si les deux
    divergent, c'est soit une cote qui a change, soit un chiffre de
    l'annexe qui n'a pas suivi -- dans les deux cas mu sortira faux et
    personne ne le verra sur le resultat seul.
    """
    problems = []
    dk = abs(cfg.MU_K - cfg.MU_K_DERIVED)
    da = abs(cfg.MU_ALPHA_DEG - cfg.MU_ALPHA_DEG_DERIVED)
    if dk > 0.01:
        problems.append(
            "k = {:.4f} (annexe A) contre {:.4f} recalcule depuis Rb = {:.1f}, "
            "rc = {:.1f}, h = {:.1f} mm".format(
                cfg.MU_K, cfg.MU_K_DERIVED, 1e3 * cfg.BALL_RADIUS_M,
                1e3 * cfg.ORING_CORD_R_M, 1e3 * cfg.ORING_SPACING_M))
    if da > 0.5:
        problems.append(
            "alpha = {:.2f} deg (annexe A) contre {:.2f} deg recalcule"
            .format(cfg.MU_ALPHA_DEG, cfg.MU_ALPHA_DEG_DERIVED))
    return problems


def check_config():
    """Coherence des seuils entre eux, avant tout mouvement."""
    problems = check_geometry_consistency()

    n_win = int(round(cfg.SETTLE_WINDOW_S * cfg.LOOP_HZ))
    if n_win < 8:
        problems.append(
            "SETTLE_WINDOW_S = {:.2f} s ne fait que {} images a {:.0f} Hz : "
            "la pente y est trop bruitee pour un seuil a {:.1f} deg/s"
            .format(cfg.SETTLE_WINDOW_S, n_win, cfg.LOOP_HZ,
                    cfg.SETTLE_RATE_DEG_S))

    # Bruit theorique de la pente sur la fenetre, a comparer au seuil.
    sigma_slope = (cfg.SIGMA_PSI_DEG
                   * math.sqrt(12.0 / (n_win * (n_win ** 2 - 1)))
                   * cfg.LOOP_HZ)
    if cfg.SETTLE_RATE_DEG_S < 3.0 * sigma_slope:
        problems.append(
            "SETTLE_RATE_DEG_S = {:.1f} deg/s est a {:.1f} sigma du bruit de "
            "pente ({:.2f} deg/s) : la stabilisation ne sera jamais declaree"
            .format(cfg.SETTLE_RATE_DEG_S,
                    cfg.SETTLE_RATE_DEG_S / sigma_slope, sigma_slope))

    if cfg.SLIP_RESID_DEG < 3.0 * cfg.SIGMA_PSI_DEG:
        problems.append(
            "SLIP_RESID_DEG = {:.1f} deg n'est qu'a {:.1f} sigma du bruit de "
            "psi ({:.2f} deg) : faux decrochages garantis"
            .format(cfg.SLIP_RESID_DEG,
                    cfg.SLIP_RESID_DEG / cfg.SIGMA_PSI_DEG, cfg.SIGMA_PSI_DEG))

    if cfg.SLIP_RESID_DEG > 2.0 * cfg.STEP_DEG:
        problems.append(
            "SLIP_RESID_DEG = {:.1f} deg vaut plus de deux pas ({:.1f} deg) : "
            "le decrochage serait date au moins deux pas trop tard"
            .format(cfg.SLIP_RESID_DEG, cfg.STEP_DEG))

    # L'acceleration de la rampe entre dans l'equilibre de la bille au meme
    # titre que la gravite : la manip n'est quasi-statique que si elle est
    # negligeable devant la limite de non-glissement.
    u_slip_bottom = cfg.U_SLIP_MAX
    if cfg.TRAP_ACCEL_RAD_S2 > 0.05 * u_slip_bottom:
        problems.append(
            "TRAP_ACCEL_RAD_S2 = {:.2f} rad/s^2 fait {:.0f} % de la limite de "
            "non-glissement au fond ({:.0f}) : la mesure n'est plus "
            "quasi-statique".format(cfg.TRAP_ACCEL_RAD_S2,
                                    100.0 * cfg.TRAP_ACCEL_RAD_S2 / u_slip_bottom,
                                    u_slip_bottom))

    return problems


# ============================================================================
# ODRIVE EN COMMANDE DE POSITION
# ============================================================================

def connect_position_mode():
    """Connexion, verifications, puis mode POSITION + TRAP_TRAJ.

    T1/T2/T3 arment en vitesse (bc.connect_and_prepare) ; E4 a besoin de
    position, d'ou une fonction propre plutot qu'un drapeau de plus dans la
    fonction commune -- la seule chose que les deux partagent est
    preflight_check, qui est reutilisee telle quelle.
    """
    print("[odrv] connexion ...")
    odrv0 = odrive.find_any()
    ax = odrv0.axis0
    print("[odrv] connecte (serial {}, fw {}.{}.{})".format(
        odrv0.serial_number, odrv0.fw_version_major,
        odrv0.fw_version_minor, odrv0.fw_version_revision))

    bc.preflight_check(odrv0, ax, cfg)

    i_soft = float(ax.config.motor.current_soft_max)
    if i_soft > cfg.I_MAX_A + 1e-6:
        raise SystemExit(
            "limite de courant ODrive a {:.1f} A, au-dela des {:.1f} A des "
            "campagnes precedentes\n(TAU_MAX / KT). E4 ne reconfigure pas le "
            "driver : la ramener a {:.1f} A dans odrivetool, sinon cet essai "
            "n'est pas comparable aux autres."
            .format(i_soft, cfg.I_MAX_A, cfg.I_MAX_A))
    print("[odrv] courant limite : {:.1f} A (plafond campagne {:.1f} A)"
          .format(i_soft, cfg.I_MAX_A))

    if ax.current_state != AXIS_STATE_CLOSED_LOOP_CONTROL:
        print("[odrv] calibration de l'offset encodeur ...")
        ax.requested_state = AXIS_STATE_ENCODER_OFFSET_CALIBRATION
        while ax.current_state != AXIS_STATE_IDLE:
            time.sleep(0.1)
        if ax.active_errors != 0 or ax.disarm_reason != 0:
            raise RuntimeError(
                "calibration echouee : active=0x{:X}, disarm=0x{:X}"
                .format(ax.active_errors, ax.disarm_reason))

    vel_tps = cfg.TRAP_VEL_RAD_S / (2.0 * math.pi)
    acc_tps2 = cfg.TRAP_ACCEL_RAD_S2 / (2.0 * math.pi)

    # L'ORDRE COMPTE, ET IL N'EST PAS ANODIN. Si l'axe est deja ARME en
    # mode vitesse -- ce qui est le cas apres une serie T1/T2 lancee avec
    # KEEP_ARMED=1 -- basculer control_mode en POSITION fait prendre effet
    # a l'input_pos courant, qui vaut ce qu'une session precedente y a
    # laisse, typiquement 0. Le cerceau partirait rejoindre cette position
    # avant meme que la consigne d'origine soit ecrite. On fige donc
    # input_pos sur la position REELLE d'abord ; en mode vitesse cette
    # ecriture est sans effet, et elle est deja la bonne au basculement.
    theta0_turns = float(ax.pos_estimate)
    ax.controller.input_pos = theta0_turns

    ax.trap_traj.config.vel_limit = vel_tps
    ax.trap_traj.config.accel_limit = acc_tps2
    ax.trap_traj.config.decel_limit = acc_tps2
    ax.controller.config.vel_limit = max(
        float(ax.controller.config.vel_limit), 3.0 * vel_tps)

    ax.controller.config.input_mode = INPUT_MODE_TRAP_TRAJ
    ax.controller.config.control_mode = CONTROL_MODE_POSITION_CONTROL
    ax.controller.input_pos = theta0_turns

    if ax.current_state != AXIS_STATE_CLOSED_LOOP_CONTROL:
        ax.requested_state = AXIS_STATE_CLOSED_LOOP_CONTROL
        time.sleep(0.5)

    if ax.active_errors != 0 or ax.disarm_reason != 0:
        raise RuntimeError("ODrive en erreur apres armement : "
                           "active=0x{:X}, disarm=0x{:X}"
                           .format(ax.active_errors, ax.disarm_reason))

    print("[odrv] mode POSITION / TRAP_TRAJ : {:.3f} rad/s, {:.2f} rad/s^2"
          .format(cfg.TRAP_VEL_RAD_S, cfg.TRAP_ACCEL_RAD_S2))
    print("[odrv] origine theta0 = {:.5f} tr".format(theta0_turns))
    return odrv0, ax, theta0_turns


def check_errors(ax, where):
    a, d = ax.active_errors, ax.disarm_reason
    if a != 0 or d != 0:
        raise RuntimeError("ODrive en erreur {} : active=0x{:X}, disarm=0x{:X}"
                           .format(where, a, d))


# ============================================================================
# FLUX DE MESURE
# ============================================================================

class Trace:
    """Journal continu de l'essai : une ligne par image."""

    def __init__(self, n_max):
        self.t = np.full(n_max, np.nan)
        self.theta = np.full(n_max, np.nan)      # encodeur [rad]
        self.theta_cmd = np.full(n_max, np.nan)  # consigne  [rad]
        self.psi = np.full(n_max, np.nan)        # camera, offset retranche
        self.ok = np.zeros(n_max, dtype=bool)
        self.step_idx = np.full(n_max, -1, dtype=np.int32)
        self.phase = np.zeros(n_max, dtype=np.int8)   # 0 repos 1 attente 2 enreg
        self.n = 0

    def append(self, t, theta, theta_cmd, psi, ok, step_idx, phase):
        i = self.n
        if i >= self.t.size:            # journal plein : on ne perd que la fin
            return
        self.t[i] = t
        self.theta[i] = theta
        self.theta_cmd[i] = theta_cmd
        self.psi[i] = psi if ok else np.nan
        self.ok[i] = ok
        self.step_idx[i] = step_idx
        self.phase[i] = phase
        self.n = i + 1

    def arrays(self):
        s = slice(0, self.n)
        return {"t": self.t[s], "theta": self.theta[s],
                "theta_cmd": self.theta_cmd[s], "psi": self.psi[s],
                "ok": self.ok[s], "step_idx": self.step_idx[s],
                "phase": self.phase[s]}


def _slope_deg_s(t, psi):
    """Pente de psi par MOINDRES CARRES sur la fenetre, en deg/s.

    Pas une difference arriere : a sigma_psi = 0.33 deg et dt = 20 ms, la
    difference arriere donne 23 deg/s de bruit sur une bille immobile,
    contre 0.64 deg/s pour la pente sur 20 points.
    """
    if t.size < 4:
        return float("inf")
    tc = t - t.mean()
    den = float(np.dot(tc, tc))
    if den <= 0.0:
        return float("inf")
    return math.degrees(float(np.dot(tc, psi - psi.mean()) / den))


class Sampler:
    """Une image -> (t, theta, psi) journalise, avec fenetre glissante.

    Concentre les trois acces materiels d'un tour (camera, pos_estimate,
    horloge) pour que la cadence soit celle de la camera et pas celle du
    code appelant.
    """

    def __init__(self, ax, picam2, det, roi, centre, theta0_turns, trace, t0):
        self.ax, self.picam2, self.det = ax, picam2, det
        self.roi, self.centre = roi, centre
        self.theta0 = theta0_turns
        self.trace, self.t0 = trace, t0
        self.psi_offset = 0.0
        self.lost = 0
        self.win_t = []
        self.win_psi = []

    def tick(self, step_idx, phase, theta_cmd):
        psi_raw = bc._read_psi_raw(self.picam2, self.det, self.roi, self.centre)
        t = time.perf_counter() - self.t0
        theta = (float(self.ax.pos_estimate) - self.theta0) * 2.0 * math.pi

        ok = psi_raw is not None
        psi = (psi_raw - self.psi_offset) if ok else float("nan")
        self.lost = 0 if ok else self.lost + 1
        if self.lost > cfg.ABORT_LOST_FRAMES:
            raise RuntimeError(
                "bille perdue sur {} images consecutives ({:.1f} s). "
                "Detection ou eclairage.".format(self.lost,
                                                 self.lost / cfg.LOOP_HZ))

        self.trace.append(t, theta, theta_cmd, psi, ok, step_idx, phase)
        if ok:
            self.win_t.append(t)
            self.win_psi.append(psi)
            while self.win_t and (t - self.win_t[0]) > cfg.SETTLE_WINDOW_S:
                self.win_t.pop(0)
                self.win_psi.pop(0)
        return t, theta, psi, ok

    def rate_deg_s(self):
        return _slope_deg_s(np.asarray(self.win_t), np.asarray(self.win_psi))

    def clear_window(self):
        self.win_t.clear()
        self.win_psi.clear()


def wait_settled(sampler, step_idx, theta_cmd):
    """Attend |dpsi/dt| < SETTLE_RATE_DEG_S pendant SETTLE_HOLD_S.

    Retourne (settled, t_settle_s). Un pas non stabilise n'est PAS jete :
    il est enregistre avec settled = False, et l'analyse decide. Jeter une
    mesure parce qu'elle est genante est la facon de fabriquer un mu.
    """
    t_start = time.perf_counter()
    t_ok_since = None
    sampler.clear_window()

    while True:
        sampler.tick(step_idx, 1, theta_cmd)
        elapsed = time.perf_counter() - t_start

        rate = sampler.rate_deg_s()
        if abs(rate) < cfg.SETTLE_RATE_DEG_S and len(sampler.win_t) >= 8:
            if t_ok_since is None:
                t_ok_since = time.perf_counter()
            elif time.perf_counter() - t_ok_since >= cfg.SETTLE_HOLD_S:
                return True, elapsed
        else:
            t_ok_since = None

        if elapsed > cfg.SETTLE_TIMEOUT_S:
            return False, elapsed


def record_step(sampler, step_idx, theta_cmd):
    """Enregistre (theta, psi) sur RECORD_S. Retourne les agregats."""
    th, ps = [], []
    t_start = time.perf_counter()
    while time.perf_counter() - t_start < cfg.RECORD_S:
        _, theta, psi, ok = sampler.tick(step_idx, 2, theta_cmd)
        th.append(theta)
        if ok:
            ps.append(psi)

    if len(ps) < 5:
        raise RuntimeError(
            "seulement {} detections sur les {:.0f} ms d'enregistrement du "
            "pas {}".format(len(ps), 1e3 * cfg.RECORD_S, step_idx))

    return {
        "theta": float(np.mean(th)),
        "theta_std": float(np.std(th)),
        "psi": float(np.mean(ps)),
        "psi_std": float(np.std(ps)),
        "n": len(ps),
    }


# ============================================================================
# DETECTION DU DECROCHAGE
# ============================================================================

def fit_follow_line(theta, psi):
    """psi = a*theta + b sur les premiers pas. Retourne (a, b).

    L'ajustement absorbe l'offset residuel de la camera et le sens de
    comptage de l'encodeur, qui n'a aucune raison de coincider avec le sens
    trigonometrique de l'image. C'est aussi le controle de sanite du
    montage : a hors de [FIT_SLOPE_MIN, FIT_SLOPE_MAX] veut dire que la
    bille ne suit pas le cerceau, et rien de ce qui suit n'a de sens.
    """
    A = np.column_stack([np.asarray(theta, float), np.ones(len(theta))])
    coef, *_ = np.linalg.lstsq(A, np.asarray(psi, float), rcond=None)
    return float(coef[0]), float(coef[1])


def lag_deg(psi, theta, a, b, direction):
    """Retard de la bille sur le cerceau, en degres, positif = elle decroche.

    Le residu psi - (a*theta + b) est signe : la bille qui decroche reste
    EN ARRIERE du cerceau, donc son residu est de signe oppose au sens de
    progression de psi, lui-meme donne par sign(a) * direction.
    """
    resid = psi - (a * theta + b)
    return -math.degrees(resid) * math.copysign(1.0, a) * direction


# ============================================================================
# UN ESSAI
# ============================================================================

def measure_psi_offset(sampler, duration):
    """Offset camera, bille au repos au fond : la physique impose psi = 0.

    Memes criteres d'acceptation que bc.arm() etape 1, et pour la meme
    raison : un offset de plus de ARM_PSI_MAX_DEG veut dire que la bille
    n'est pas au fond ou que HOOP_CENTRE_PX est faux, pas qu'il faut le
    retrancher.
    """
    t_list, psi_list = [], []
    t_start = time.perf_counter()
    sampler.psi_offset = 0.0
    while time.perf_counter() - t_start < duration:
        t, _, psi, ok = sampler.tick(-1, 0, 0.0)
        if ok:
            t_list.append(t)
            psi_list.append(psi)

    if len(psi_list) < 20:
        raise RuntimeError(
            "bille detectee sur {} images seulement pendant la mesure "
            "d'offset".format(len(psi_list)))

    st = bc.psi_statistics(np.asarray(t_list), np.asarray(psi_list), cfg.F_N_HZ)
    print("      offset camera      : {:+6.2f} deg  (limite {:.1f})"
          .format(math.degrees(st["offset"]), cfg.ARM_PSI_MAX_DEG))
    print("      oscillation {:.2f} Hz : {:6.2f} deg  (limite {:.1f})"
          .format(cfg.F_N_HZ, math.degrees(st["osc_amp"]), cfg.ARM_OSC_MAX_DEG))
    print("      bruit de detection : {:6.2f} deg  (limite {:.1f})"
          .format(math.degrees(st["sigma"]), cfg.ARM_NOISE_MAX_DEG))

    if math.degrees(st["osc_amp"]) > cfg.ARM_OSC_MAX_DEG:
        raise RuntimeError(
            "la bille oscille encore ({:.2f} deg a {:.2f} Hz) : laisser "
            "s'immobiliser".format(math.degrees(st["osc_amp"]), cfg.F_N_HZ))
    if math.degrees(st["sigma"]) > cfg.ARM_NOISE_MAX_DEG:
        raise RuntimeError(
            "detection trop bruitee ({:.2f} deg) : c'est la vision, pas la "
            "bille".format(math.degrees(st["sigma"])))
    if abs(math.degrees(st["offset"])) > cfg.ARM_PSI_MAX_DEG:
        raise RuntimeError(
            "offset camera de {:+.1f} deg : bille pas au fond, ou "
            "HOOP_CENTRE_PX faux".format(math.degrees(st["offset"])))

    sampler.psi_offset = st["offset"]
    return st


def run_trial(ax, picam2, det, roi, centre, theta0_turns, direction, trial):
    """Un essai complet dans un sens. Retourne (trace_arrays, steps, meta)."""
    print("\n" + "-" * 70)
    print("  essai {} -- sens {:+d}".format(trial, direction))
    print("-" * 70)

    # --- retour a l'origine et immobilisation ---------------------------
    ax.controller.input_pos = theta0_turns
    time.sleep(cfg.REZERO_SETTLE_S)
    check_errors(ax, "au retour a l'origine")

    n_max = int(cfg.LOOP_HZ * (cfg.ARM_SETTLE_S + 5.0
                               + (cfg.ABORT_THETA_DEG / cfg.STEP_DEG)
                               * (cfg.SETTLE_TIMEOUT_S + cfg.RECORD_S)) + 500)
    trace = Trace(n_max)
    t0 = time.perf_counter()
    sampler = Sampler(ax, picam2, det, roi, centre, theta0_turns, trace, t0)

    print("[e4] mesure de l'offset camera ({:.1f} s) ...".format(cfg.ARM_SETTLE_S))
    arm_stats = measure_psi_offset(sampler, cfg.ARM_SETTLE_S)

    # --- pas 0 : etat de depart, ancre de l'ajustement -------------------
    steps = [record_step(sampler, 0, 0.0)]
    steps[0]["settled"] = True
    steps[0]["theta_cmd"] = 0.0
    steps[0]["t_settle_s"] = 0.0

    theta_cmd = 0.0
    a = b = float("nan")
    slip_run = 0
    lag_prev = -float("inf")
    theta_star = psi_star = float("nan")
    step_star = -1
    extra_left = None
    stop_reason = "termine"

    k = 0
    while True:
        k += 1
        theta_cmd = direction * cfg.STEP_DEG * k
        if abs(theta_cmd) > cfg.ABORT_THETA_DEG:
            stop_reason = "borne theta atteinte sans decrochage"
            break

        ax.controller.input_pos = theta0_turns + math.radians(theta_cmd) / (2.0 * math.pi)
        settled, t_settle = wait_settled(sampler, k, math.radians(theta_cmd))
        rec = record_step(sampler, k, math.radians(theta_cmd))
        rec["settled"] = bool(settled)
        rec["theta_cmd"] = math.radians(theta_cmd)
        rec["t_settle_s"] = float(t_settle)
        steps.append(rec)
        check_errors(ax, "pendant l'increment")

        # --- ajustement de la droite de suivi, une fois pour l'essai -----
        if k == cfg.N_FIT_STEPS:
            a, b = fit_follow_line([s["theta"] for s in steps],
                                   [s["psi"] for s in steps])
            print("[e4] droite de suivi : psi = {:+.3f} theta {:+.3f} rad"
                  .format(a, b))
            if not (cfg.FIT_SLOPE_MIN <= abs(a) <= cfg.FIT_SLOPE_MAX):
                stop_reason = "pente de suivi aberrante"
                raise RuntimeError(
                    "pente psi/theta = {:+.3f}, hors de [{:.1f}, {:.1f}] : la "
                    "bille ne suit pas le cerceau sur les {} premiers pas. "
                    "Verifier qu'elle est bien dans la gorge et que le sens "
                    "de l'encodeur est celui attendu."
                    .format(a, cfg.FIT_SLOPE_MIN, cfg.FIT_SLOPE_MAX,
                            cfg.N_FIT_STEPS))

        # --- surveillance du residu -------------------------------------
        if k >= cfg.N_FIT_STEPS:
            lag = lag_deg(rec["psi"], rec["theta"], a, b, direction)
            rec["lag_deg"] = lag
            print("  pas {:3d} : theta = {:+6.2f} deg  psi = {:+6.2f} deg  "
                  "retard = {:+5.2f} deg  {}{}"
                  .format(k, math.degrees(rec["theta"]),
                          math.degrees(rec["psi"]), lag,
                          "stab" if settled else "NON STAB",
                          "" if extra_left is None else "  (apres decrochage)"))

            if extra_left is None:
                if lag > cfg.SLIP_RESID_DEG and lag >= lag_prev:
                    slip_run += 1
                else:
                    slip_run = 0
                lag_prev = lag

                if slip_run >= cfg.SLIP_CONFIRM_STEPS:
                    # theta* = dernier pas ou la bille suivait encore, donc
                    # celui d'avant le PREMIER des pas confirmes.
                    step_star = k - cfg.SLIP_CONFIRM_STEPS
                    theta_star = steps[step_star]["theta"]
                    psi_star = steps[step_star]["psi"]
                    print("\n[e4] DECROCHAGE confirme au pas {} : theta* = "
                          "{:+.2f} deg, psi* = {:+.2f} deg -> mu = {:.3f}"
                          .format(k, math.degrees(theta_star),
                                  math.degrees(psi_star),
                                  cfg.mu_required(psi_star)))
                    extra_left = cfg.EXTRA_STEPS_AFTER_SLIP
                    stop_reason = "decrochage detecte"
            else:
                extra_left -= 1
                if extra_left <= 0:
                    break
        else:
            print("  pas {:3d} : theta = {:+6.2f} deg  psi = {:+6.2f} deg  "
                  "(ajustement){}".format(k, math.degrees(rec["theta"]),
                                          math.degrees(rec["psi"]),
                                          "" if settled else "  NON STAB"))

    ax.controller.input_pos = theta0_turns

    meta = {
        "direction": int(direction),
        "theta_star_rad": float(theta_star),
        "theta_star_deg": math.degrees(theta_star) if np.isfinite(theta_star) else None,
        "psi_star_rad": float(psi_star),
        "psi_star_deg": math.degrees(psi_star) if np.isfinite(psi_star) else None,
        "mu_req": cfg.mu_required(psi_star) if np.isfinite(psi_star) else None,
        "step_star": int(step_star),
        "n_steps": len(steps) - 1,
        "fit_slope": a, "fit_offset_rad": b,
        "stop_reason": stop_reason,
        "slip_detected": bool(np.isfinite(theta_star)),
        "arm_offset_deg": math.degrees(arm_stats["offset"]),
        "arm_noise_deg": math.degrees(arm_stats["sigma"]),
        "psi_offset_rad": float(sampler.psi_offset),
    }
    return trace.arrays(), steps, meta


# ============================================================================
# SAUVEGARDE
# ============================================================================

def save_trial(outdir, trial, direction, trace, steps, trial_meta, roi, centre):
    os.makedirs(outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(outdir, "E4_{}_mu_{}_t{:02d}_{}.npz".format(
        cfg.HOOP, "pos" if direction > 0 else "neg", trial, stamp))

    def col(key, dtype=float):
        return np.asarray([s.get(key, np.nan) for s in steps], dtype=dtype)

    meta = {
        "experiment": "E4",
        "version": 3,
        "timestamp": stamp,
        "bench_version": 3,
        "hoop": cfg.HOOP,
        "trial": trial,
        # --- protocole ---
        "step_deg": cfg.STEP_DEG,
        "trap_vel_rad_s": cfg.TRAP_VEL_RAD_S,
        "trap_accel_rad_s2": cfg.TRAP_ACCEL_RAD_S2,
        "settle_window_s": cfg.SETTLE_WINDOW_S,
        "settle_rate_deg_s": cfg.SETTLE_RATE_DEG_S,
        "settle_hold_s": cfg.SETTLE_HOLD_S,
        "settle_timeout_s": cfg.SETTLE_TIMEOUT_S,
        "record_s": cfg.RECORD_S,
        "slip_resid_deg": cfg.SLIP_RESID_DEG,
        "slip_confirm_steps": cfg.SLIP_CONFIRM_STEPS,
        "n_fit_steps": cfg.N_FIT_STEPS,
        "abort_theta_deg": cfg.ABORT_THETA_DEG,
        # --- modele du contact ---
        "mu_k": cfg.MU_K,
        "mu_alpha_deg": cfg.MU_ALPHA_DEG,
        "mu_k_derived": cfg.MU_K_DERIVED,
        "mu_alpha_deg_derived": cfg.MU_ALPHA_DEG_DERIVED,
        "mu_formula": "mu = ((k-1)/k)*cos(alpha)*tan(psi)  [annexe A]",
        # --- banc, identique aux autres campagnes ---
        "tau_max": cfg.TAU_MAX, "kt": cfg.KT, "i_max_a": cfg.I_MAX_A,
        "loop_hz": cfg.LOOP_HZ,
        "cam_latency_s": cfg.CAM_LATENCY_S,
        "sigma_psi_deg": cfg.SIGMA_PSI_DEG,
        "hoop_centre_px": list(cfg.HOOP_CENTRE_PX),
        "hoop_radius_px": cfg.HOOP_RADIUS_PX,
        "roi": [int(v) for v in roi],
        "centre_px": [float(v) for v in centre],
        "camera": {"size": list(cfg.SIZE), "fps": cfg.FPS,
                   "exposure_us": cfg.EXPOSURE_US, "gain": cfg.GAIN_CAM},
        "detector": {"threshold": cfg.THRESHOLD, "downsample": cfg.DOWNSAMPLE,
                     "tracking_window": cfg.TRACKING_WIN,
                     "ball_size": list(cfg.BALL_SIZE)},
    }
    meta.update(trial_meta)

    np.savez_compressed(
        path,
        meta=json.dumps(meta),
        # --- flux continu ---
        t=trace["t"], theta=trace["theta"], theta_cmd=trace["theta_cmd"],
        psi=trace["psi"], ok=trace["ok"], step_idx=trace["step_idx"],
        phase=trace["phase"],
        # --- agregats par pas ---
        step_theta=col("theta"), step_theta_std=col("theta_std"),
        step_psi=col("psi"), step_psi_std=col("psi_std"),
        step_theta_cmd=col("theta_cmd"), step_lag_deg=col("lag_deg"),
        step_settled=np.asarray([bool(s["settled"]) for s in steps]),
        step_n=col("n"), step_t_settle_s=col("t_settle_s"),
    )
    print("[out] {}".format(path))
    return path


# ============================================================================
# MAIN
# ============================================================================

def main():
    p = argparse.ArgumentParser(
        description="E4 : mesure de mu par montee quasi-statique du cerceau")
    p.add_argument("--outdir", default=cfg.OUTDIR)
    p.add_argument("--trials", type=int, default=cfg.N_TRIALS_PER_DIR,
                   help="essais par sens (defaut {})".format(cfg.N_TRIALS_PER_DIR))
    p.add_argument("--direction", type=int, choices=[1, -1], default=None,
                   help="ne faire qu'un sens")
    p.add_argument("--resume", type=int, default=1,
                   help="numero du premier essai (reprise de campagne)")
    p.add_argument("--no-prompt", action="store_true",
                   help="enchaine les essais sans attendre l'operateur")
    p.add_argument("--dry-run", action="store_true",
                   help="verifie la configuration, ne touche a rien")
    args = p.parse_args()

    directions = (args.direction,) if args.direction else cfg.DIRECTIONS

    print("=" * 70)
    print("  E4 -- coefficient de frottement statique bille / O-rings")
    print("=" * 70)
    print("  cerceau         : {}".format(cfg.HOOP))
    print("  contact         : k = {:.4f}, alpha = {:.2f} deg"
          .format(cfg.MU_K, cfg.MU_ALPHA_DEG))
    print("                    (recalcules : {:.4f}, {:.2f} deg)"
          .format(cfg.MU_K_DERIVED, cfg.MU_ALPHA_DEG_DERIVED))
    print("  pas             : {:.2f} deg a {:.2f} rad/s, {:.2f} rad/s^2"
          .format(cfg.STEP_DEG, cfg.TRAP_VEL_RAD_S, cfg.TRAP_ACCEL_RAD_S2))
    print("  stabilisation   : |dpsi/dt| < {:.1f} deg/s pendant {:.2f} s "
          "(fenetre {:.2f} s)".format(cfg.SETTLE_RATE_DEG_S, cfg.SETTLE_HOLD_S,
                                      cfg.SETTLE_WINDOW_S))
    print("  decrochage      : retard > {:.1f} deg sur {} pas ({:.1f} sigma "
          "du bruit psi)".format(cfg.SLIP_RESID_DEG, cfg.SLIP_CONFIRM_STEPS,
                                 cfg.SLIP_RESID_DEG / cfg.SIGMA_PSI_DEG))
    print("  borne theta     : {:.0f} deg".format(cfg.ABORT_THETA_DEG))
    for mu in (0.3, 0.5, 1.0):
        print("  si mu = {:.1f}      : decrochage attendu a psi* = {:.1f} deg"
              .format(mu, math.degrees(cfg.psi_star_for_mu(mu))))
    print("  mu mesurable    : jusqu'a {:.2f} (au-dela, borne theta atteinte)"
          .format(cfg.mu_required(math.radians(cfg.ABORT_THETA_DEG - cfg.STEP_DEG))))
    print("  essais          : {} x {} sens".format(args.trials, len(directions)))

    problems = check_config()
    if problems:
        print("\n[!!] configuration incoherente :")
        for q in problems:
            print("     - {}".format(q))
        raise SystemExit("corriger e4_config.py")
    print("\n  configuration coherente")

    if args.dry_run:
        print("  --dry-run : rien n'a ete envoye au materiel")
        return

    picam2 = det = ax = None
    results = []
    try:
        picam2, det, roi, centre = bc.setup_camera(cfg)
        odrv0, ax, theta0_turns = connect_position_mode()

        trial = args.resume
        for direction in directions:
            for _ in range(args.trials):
                if not args.no_prompt:
                    input("\n[e4] essai {} sens {:+d} -- poser la bille au "
                          "fond, ENTREE pour lancer ...".format(trial, direction))

                # UN ESSAI RATE NE DOIT PAS TUER LA CAMPAGNE. Les echecs
                # attendus -- bille encore en oscillation au demarrage,
                # detection perdue sur quelques images, pente de suivi
                # aberrante -- sont des conditions d'essai, pas des pannes
                # du banc. En mode --no-prompt il n'y a personne au clavier
                # pour relancer, et perdre les huit essais suivants parce
                # que le troisieme est parti trop tot serait absurde. Un
                # RuntimeError donne donc UN nouvel essai, puis on passe au
                # suivant. Ce qui casse vraiment (ODrive en erreur,
                # deconnexion) remonte par les autres types d'exception et
                # arrete tout, comme il se doit.
                for attempt in (1, 2):
                    try:
                        trace, steps, tmeta = run_trial(
                            ax, picam2, det, roi, centre, theta0_turns,
                            direction, trial)
                        save_trial(args.outdir, trial, direction, trace, steps,
                                   tmeta, roi, centre)
                        results.append(tmeta)
                        break
                    except RuntimeError as exc:
                        print("\n[e4] essai {} interrompu : {}".format(trial, exc))
                        if attempt == 1:
                            print("[e4] nouvelle tentative du meme essai apres "
                                  "{:.0f} s d'immobilisation ..."
                                  .format(cfg.REZERO_SETTLE_S))
                        else:
                            print("[e4] essai {} ABANDONNE, on passe au suivant."
                                  .format(trial))
                trial += 1

    finally:
        if ax is not None:
            try:
                ax.controller.input_pos = float(ax.pos_estimate)
                time.sleep(0.3)
                ax.requested_state = AXIS_STATE_IDLE
            except Exception:
                pass
        try:
            picam2.stop()
            picam2.close()
        except Exception:
            pass

    # --- resume de campagne ---------------------------------------------
    mus = [r["mu_req"] for r in results if r["mu_req"] is not None]
    print("\n" + "=" * 70)
    print("  {} essais, {} avec decrochage detecte".format(len(results), len(mus)))
    if mus:
        print("  mu = {:.3f} +/- {:.3f}  (moyenne +/- ecart-type)"
              .format(float(np.mean(mus)), float(np.std(mus, ddof=1))
                      if len(mus) > 1 else 0.0))
        print("  litterature / hypothese du modele : 0.50")
    print("  analyse detaillee : python3 e4_analyze.py {}".format(args.outdir))


if __name__ == "__main__":
    main()
