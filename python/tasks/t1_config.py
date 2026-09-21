"""
t1_config_v3.py -- Réglages de la stabilisation T1.

Le hoop est choisi avec la variable d'environnement :

    T1_HOOP=outer
    T1_HOOP=inner

La calibration encodeur est contrôlée avec :

    CALIBRATE_ENCODER=1

Par défaut, aucune calibration automatique n'est demandée.
"""

import math
import os

import numpy as np


# ============================================================================
# CONFIGURATION DU HOOP
# ============================================================================

HOOP = os.environ.get("T1_HOOP", "outer").lower()

if HOOP not in ("outer", "inner"):
    raise ValueError(
        "T1_HOOP doit valoir 'outer' ou 'inner', "
        f"pas {HOOP!r}"
    )


# Calibration explicite uniquement si demandée.
#
# Première calibration après mise sous tension :
#     CALIBRATE_ENCODER=1 python3 t1_run_v3.py ...
#
# Essais suivants :
#     python3 t1_run_v3.py ...
CALIBRATE_ENCODER = (os.environ.get("CALIBRATE_ENCODER", "0") == "1")
# Laisse l'axe arme entre deux essais : l'offset encodeur survit et la
# calibration n'est refaite qu'a la premiere execution de la session.
# Le cerceau reste tenu a vitesse nulle entre les essais.
KEEP_ARMED = os.environ.get("KEEP_ARMED", "0") == "1"

# ============================================================================
# COMMANDE
# ============================================================================

GAINS_OUTER = {
    "nominal": np.array([
        1.0066,
        5.9102,
        -79.8608,
        -19.9499,
    ]),

    "fallback": np.array([
        0.1007,
        1.6394,
        -12.2436,
        -11.6190,
    ]),

    "A": np.array([
        0.8,
        3.0,
        -20.0,
        -6.0,
    ]),

    "B": np.array([
        1.0,
        4.0,
        -18.0,
        -5.5,
    ]),

    "C": np.array([
        1.0,
        6.0,
        -30.0,
        -15.0,
    ]),
}


GAINS_INNER = {
    # Jeu fourni par la simulation rolling_in_inside,
    # avec inversion des signes de K3 et K4 pour le banc.
    "nominal": np.array([
        1.0066,
        6.1964,
        -87.8026,
        -19.8888,
    ]),
    
    "fallback": np.array([
        0.1007,
        1.6424,
        -2.1533,
        -9.7069,
    ]),

    # Valeurs provisoires uniquement.
    # Ne pas les considérer comme validées expérimentalement.
    "A": np.array([
        0.8,
        3.0,
        -20.0,
        -6.0,
    ]),

    "B": np.array([
        1.0,
        4.0,
        -18.0,
        -5.5,
    ]),

    "C": np.array([
        1.0,
        6.0,
        -30.0,
        -15.0,
    ]),
}


GAINS = (
    GAINS_INNER
    if HOOP == "inner"
    else GAINS_OUTER
)

DEFAULT_GAIN = "nominal"

U_EQ = 0.0

GAIN_START = 1.0
GAIN_RAMP_S = 0.0


# ============================================================================
# SOURCE DE THETA ET THETA_DOT
# ============================================================================

USE_MEASURED_THETA = False


# ============================================================================
# ACTIONNEUR
# ============================================================================

TAU_MAX = 0.5425
B_MOTEUR = 6.19e-5
I_TOTAL = 1.1193e-3
KT = 0.031

U_SLIP_MAX = 100.0
THETA_DOT_MAX = 20.0

# Constante de temps equivalente de la boucle de vitesse ODrive, vue depuis
# le Pi (variateur + moteur + retard USB). 0.0 = hypothese de boucle interne
# parfaite, c'est-a-dire le comportement anterieur au bit pres. Mesuree par
# vel_loop_id_v3.py. Laissee a zero ici tant que les campagnes T1/T2 ne sont
# pas refaites : leur ecart |vel_cmd - vel_enc| (0.15-0.40 rad/s RMS) reste
# dans le bruit, donc le filtre n'y changerait rien de mesurable.
TAU_VEL_S = 0.0


# ============================================================================
# MODELE LINEARISE DE LA BALLE
# ============================================================================

if HOOP == "outer":
    F_N_HZ = 1.3144
    ZETA = 0.033
    B2 = -0.402

else:
    F_N_HZ = 2.2157
    ZETA   = 0.0198
    B2     = -0.4720

OMEGA_N = 2.0 * math.pi * F_N_HZ
A21 = -OMEGA_N ** 2
A22 = -2.0 * ZETA * OMEGA_N


# ============================================================================
# CAMERA
# ============================================================================

SIZE = (820, 616)
FPS = 50.0
EXPOSURE_US = 18000
GAIN_CAM = 12.0

THRESHOLD = 110
DOWNSAMPLE = 8
TRACKING_WIN = 64
BALL_SIZE = (0, 150)

ROI_MARGIN = 0.05
CAM_LATENCY_S = 0.0193


# Centres exprimés en pixels dans l'image complète.
#
# Le centre outer est celui mesuré avec Hough.
# Le centre inner doit être vérifié avec une détection Hough dédiée.
HOOP_CENTRES_PX = {
    "outer": (384.60, 276.60),
    "inner": (384.60, 276.60),
}

# Rayons exprimés en pixels.
#
# Le rayon inner est une estimation initiale et doit être vérifié.
HOOP_RADII_PX = {
    "outer": 269.0,
    "inner": 116.0,
}

HOOP_CENTRE_PX = HOOP_CENTRES_PX[HOOP]
HOOP_RADIUS_PX = HOOP_RADII_PX[HOOP]

USE_HOUGH_CENTRE = True


# ============================================================================
# ESTIMATION DE PSI_DOT
# ============================================================================

PSIDOT_FILTER_HZ = 12.0

SIGMA_PSI_DEG = 0.4
Q_ACCEL = 9.0


# ============================================================================
# BOUCLE ET CADENCES ODRIVE
# ============================================================================

LOOP_HZ = 50.0
DT_NOMINAL = 1.0 / LOOP_HZ
DT_MAX = 0.060

COM_PERIOD = 1
ENC_PERIOD = 20
IQ_PERIOD = 0
ERR_PERIOD = 20


# ============================================================================
# SECURITES
# ============================================================================

ARM_SETTLE_S = 2.0
ARM_HOLD_S = 0.4

ARM_PSI_MAX_DEG = 8.0
ARM_OSC_MAX_DEG = 1.5
ARM_NOISE_MAX_DEG = 2.0
ARM_INIT_TOL_DEG = 10.0

ABORT_PSI_DEG = 100.0
ABORT_LOST_FRAMES = 5
ABORT_OVERRUN = 10
MAX_DURATION_S = 60.0

RAMP_DOWN_S = 0.8


# ============================================================================
# SORTIE
# ============================================================================

OUTDIR = (
    "data/t1_inner"
    if HOOP == "inner"
    else "data/t1_outer"
)