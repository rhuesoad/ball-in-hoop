#!/usr/bin/env python3
"""
e4_config.py -- E4 : mesure du coefficient de frottement statique bille /
O-rings par montee quasi-statique du cerceau.

CE QUE MESURE E4, ET POURQUOI IL FAUT LE MESURER
================================================
params.ball_track_friction_coeff vaut 0.5 dans ball_hoop_params.m avec la
mention "CONSERVATIVE, not measured", et U_SLIP_MAX cote banc en descend.
Toute la planification T3/T4 sous contrainte de non-glissement repose donc
sur un nombre emprunte a la litterature. E4 le remplace par une mesure.

PRINCIPE
--------
Le cerceau est commande EN POSITION (et non en vitesse comme T1/T2/T3). On
part de theta = 0, bille au repos au fond, et on incremente theta par pas
de STEP_DEG en laissant la bille se stabiliser a chaque pas.

Tant que le contact tient, la bille est portee par le cerceau et monte avec
lui : psi suit theta. Au-dela d'un certain angle la friction disponible ne
suffit plus, la bille glisse dans la gorge et decroche : un ecart
theta - psi apparait et croit. L'angle du dernier pas ou la bille suivait
encore est theta*, et l'annexe A du memoire en tire

    mu_req = ((k-1)/k) * cos(alpha) * tan(psi)

avec k le facteur d'inertie effective du contact deux rails et alpha le
demi-angle de la gorge. Voir MU_K et MU_ALPHA_DEG plus bas.

CE QUI EST HERITE DE T1, ET CE QUI NE L'EST PAS
-----------------------------------------------
Toute la chaine de mesure est celle de T1/T2/T3 : meme camera, meme
detecteur, meme ROI, meme encodeur, memes limites moteur. C'est la
condition pour que le sigma_psi caracterise en E2/commissioning s'applique
ici, et donc pour que le seuil de decrochage soit justifie plutot que
choisi.

Ce qui change est la COMMANDE : position au lieu de vitesse, et une rampe
trapezoidale volontairement tres lente. L'acceleration du cerceau entre
dans l'equilibre de la bille exactement comme la gravite ; si elle n'est
pas negligeable devant la limite de non-glissement, la manip mesure autre
chose que du statique. Voir TRAP_ACCEL_RAD_S2.
"""

import math
import os

# Banc, camera, detecteur, modele linearise, securites d'armement.
# L'import est en * comme dans t2_config_v3 / t3_config_v3 : une liste
# explicite se desynchronise en silence et le trou n'apparait qu'au premier
# essai reel.
from tasks.t1_config import *          # noqa: F401,F403,E402
from tasks import t1_config as _base        # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))


# ============================================================================
# GEOMETRIE DU CONTACT DEUX RAILS  (annexe A du memoire)
# ============================================================================
# Les deux constantes sont RECALCULEES ci-dessous a partir des cotes
# mesurees, et l'ecart avec les valeurs de l'annexe est verifie au
# lancement (e4_run.py, check_geometry_consistency). Deux sources de verite
# qui derivent l'une de l'autre sont la facon habituelle de se retrouver
# avec un mu faux de 10 % sans le voir.
MU_K         = 1.552        # [-]   k = 1 + (2/5)(Rb/r_roll)^2
MU_ALPHA_DEG = 31.7         # [deg] demi-angle du contact dans la gorge

# Cotes as-built (ball_hoop_params.m : mesures, pas des hypotheses).
BALL_RADIUS_M   = 12.5e-3     # Rb
ORING_CORD_R_M  = 1.5e-3      # rc
ORING_SPACING_M = 14.7e-3     # h

# Delta      = sqrt((Rb + rc)^2 - (h/2)^2)   offset radial du centre bille
# r_roll     = Rb * Delta / (Rb + rc)        rayon de roulement, < Rb
# cos(alpha) = Delta / (Rb + rc)             d'ou alpha
_DELTA_M  = math.sqrt((BALL_RADIUS_M + ORING_CORD_R_M) ** 2
                      - (ORING_SPACING_M / 2.0) ** 2)
_R_ROLL_M = BALL_RADIUS_M * _DELTA_M / (BALL_RADIUS_M + ORING_CORD_R_M)

MU_K_DERIVED         = 1.0 + 0.4 * (BALL_RADIUS_M / _R_ROLL_M) ** 2
MU_ALPHA_DEG_DERIVED = math.degrees(
    math.acos(_DELTA_M / (BALL_RADIUS_M + ORING_CORD_R_M)))


def mu_required(psi_rad):
    """mu_req = ((k-1)/k) * cos(alpha) * tan(psi).   [annexe A]

    psi en radians. La formule est impaire en psi et un essai en sens
    negatif doit donner le meme mu, d'ou la valeur absolue.
    """
    return ((MU_K - 1.0) / MU_K) * math.cos(math.radians(MU_ALPHA_DEG)) \
        * abs(math.tan(psi_rad))


def psi_star_for_mu(mu):
    """Reciproque de mu_required : angle de decrochage attendu pour un mu
    donne. Sert a dimensionner ABORT_THETA_DEG et a annoncer, avant de
    lancer, ou l'essai devrait decrocher si la litterature avait raison."""
    c = ((MU_K - 1.0) / MU_K) * math.cos(math.radians(MU_ALPHA_DEG))
    return math.atan(mu / c)


# ============================================================================
# COMMANDE EN POSITION
# ============================================================================

# Pas d'increment. Il fixe DIRECTEMENT la resolution sur theta*, donc sur
# mu : d(mu)/d(psi) = ((k-1)/k) cos(alpha) / cos^2(psi) vaut 0.020 par
# degre au voisinage de psi = 59 deg (le decrochage attendu si mu = 0.5),
# soit 0.03 sur mu pour un pas de 1.5 deg.
STEP_DEG = 1.5

# Rampe trapezoidale du deplacement d'un pas. LA VALEUR QUI COMPTE EST
# L'ACCELERATION : elle entre dans l'equilibre de la bille au meme titre
# que la gravite. La limite de non-glissement au fond vaut ~100 rad/s^2
# pour mu = 0.5 (slip_acceleration_limit.m), donc 0.5 rad/s^2 en represente
# 0.5 % : la manip reste quasi-statique par construction, ce n'est pas une
# esperance. Un pas de 1.5 deg a 0.15 rad/s dure alors ~0.20 s.
TRAP_VEL_RAD_S    = 0.15
TRAP_ACCEL_RAD_S2 = 0.50

# Nombre d'essais par sens de rotation, et sens testes.
N_TRIALS_PER_DIR = 5
DIRECTIONS = (+1, -1)


# ============================================================================
# STABILISATION ET ENREGISTREMENT A CHAQUE PAS
# ============================================================================

# psi_dot est estime par REGRESSION LINEAIRE sur une fenetre glissante, et
# non par difference arriere : a sigma_psi = 0.33 deg et dt = 20 ms, une
# difference arriere donnerait 23 deg/s de bruit sur une bille immobile,
# soit huit fois le seuil. La pente sur 20 points a un ecart-type de
# sigma * sqrt(12 / (N(N^2-1))) / dt = 0.64 deg/s, donc SETTLE_RATE_DEG_S
# = 3 deg/s est un seuil a 5 sigma.
SETTLE_WINDOW_S   = 0.40     # 20 images a 50 Hz
SETTLE_RATE_DEG_S = 3.0
SETTLE_HOLD_S     = 0.50     # duree pendant laquelle le critere doit tenir
SETTLE_TIMEOUT_S  = 3.0      # au-dela, le pas est marque non stabilise

# Fenetre d'enregistrement de (theta, psi) une fois la stabilisation
# acquise. 0.3 s = 15 images : l'incertitude sur la moyenne de psi vaut
# sigma / sqrt(15) = 0.085 deg, tres en dessous du pas de 1.5 deg.
RECORD_S = 0.30


# ============================================================================
# CRITERE DE DECROCHAGE
# ============================================================================

# L'ecart mesure n'est PAS |theta - psi| brut. Les deux angles n'ont ni la
# meme origine (offset camera) ni forcement le meme signe (sens de comptage
# de l'encodeur contre sens trigonometrique de l'image), et la gorge a une
# compliance qui donne une pente a legerement differente de 1. On ajuste
# donc psi = a*theta + b sur les N_FIT_STEPS premiers pas -- ou la bille
# suit par hypothese -- et on surveille le RESIDU par rapport a cette
# droite.
#
# a est aussi le controle de sanite du montage : hors de [0.5, 1.5] la
# bille ne suit pas le cerceau et l'essai n'a pas de sens.
N_FIT_STEPS   = 5
FIT_SLOPE_MIN = 0.5
FIT_SLOPE_MAX = 1.5

# Seuil a 6 sigma sur le bruit de psi (0.33 deg), confirme sur deux pas
# consecutifs a residu croissant : une image aberrante ne suffit pas a
# declencher. theta* est l'angle du DERNIER pas avant le premier des deux.
SLIP_RESID_DEG     = 2.0
SLIP_CONFIRM_STEPS = 2

# Pas supplementaires poursuivis apres detection, uniquement pour que le
# graphe theta-psi montre la divergence au lieu de s'arreter dessus. Ils
# n'entrent dans aucun calcul.
EXTRA_STEPS_AFTER_SLIP = 3


# ============================================================================
# SECURITES
# ============================================================================

# Borne dure sur theta. Si mu = 0.5, le decrochage est attendu vers
# psi* = atan(mu * k / ((k-1) cos alpha)) = 58.8 deg ; a 90 deg on est bien
# au-dela sans jamais mettre la bille au-dessus de l'horizontale du
# cerceau. Atteindre cette borne SANS decrochage est un resultat en soi
# (mu > 0.90) et l'essai est marque comme tel, pas jete.
ABORT_THETA_DEG = 90.0

# Courant : meme limite que les scripts existants. TAU_MAX / KT = 17.5 A.
# Le script REFUSE de partir si l'ODrive est configure au-dessus, il ne
# reconfigure rien : ecrire dans motor.config depuis un script d'essai est
# la facon sure de laisser le banc dans un etat different de celui des
# campagnes precedentes.
I_MAX_A = _base.TAU_MAX / _base.KT

# Detection : au-dela, la bille est perdue et le pas est invalide.
ABORT_LOST_FRAMES = 25          # 0.5 s sans detection

# Retour a l'origine entre deux essais : duree laissee a la bille pour
# redescendre et s'immobiliser au fond.
REZERO_SETTLE_S = 4.0


# ============================================================================
# SORTIE
# ============================================================================

OUTDIR = os.path.join(_HERE, "data", "e4_mu")
