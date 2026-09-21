"""
t2_config_v3.py -- Reglages de la poursuite T2 (cerceau exterieur).

T2 partage EXACTEMENT le meme banc, le meme modele et les memes securites
que T1 : tout est importe de t1_config_v3 pour qu'il n'existe qu'une seule
source de verite. Seuls les reglages propres a la reference variable sont
definis ici.
"""

import math

from tasks.t1_config import *          # noqa: F401,F403  (banc, modele, securites)
from tasks import t1_config as _base

# ========================================================================== #
#  REFERENCE
# ========================================================================== #

REF_TYPE    = "sine"            # "sine" | "step" | "triangle"
REF_AMP_DEG = 5.0               # amplitude [deg]
REF_FREQ_HZ = 1.0               # frequence [Hz]

REF_FADE_S = 2.0                # [s]

# ========================================================================== #
#  ENVELOPPE DE POURSUITE
# ========================================================================== #

def reference_effort(amp_deg, freq_hz):
    """Effort demande au cerceau pour une sinusoide (A, f), sans le suivi.

    Retourne (max|u_r|, max|theta_dot_r|, max|theta_r|) en unites SI.
    """
    A = math.radians(amp_deg)
    w = 2.0 * math.pi * freq_hz
    P = -A * (w ** 2 + _base.A21) / _base.B2
    Q = -A * _base.A22 * w / _base.B2
    amp = math.hypot(P, Q)
    return amp, amp / w, amp / w ** 2


def max_amplitude_deg(freq_hz):
    """Amplitude maximale de psi_r a la frequence donnee, sous la contrainte
    THETA_DOT_MAX. Lineaire en A, donc obtenue par simple mise a l'echelle."""
    _, theta_dot_1deg, _ = reference_effort(1.0, freq_hz)
    if theta_dot_1deg <= 0.0:
        return float("inf")
    return _base.THETA_DOT_MAX / theta_dot_1deg


# ========================================================================== #
#  SORTIE
# ========================================================================== #

OUTDIR = "data/t2_v3"
