#!/usr/bin/env python3
"""
t4_config_v3.py -- T4 : flying ball, pompage puis largage puis rattrapage
                        sur la face interne du cerceau interieur.

Meme principe que t3_config_v3.py : on herite de t1_config_v3 tout ce qui
decrit le BANC et on ne redefinit que ce que la manoeuvre change. Rien de
la chaine d'acquisition n'est touche.

CE QUE T4 AJOUTE PAR RAPPORT A T3
=================================
T3 est UNE trajectoire et UN gain variant. T4 est TROIS lois de commande
enchainees par un test sur l'etat mesure :

  phase 1  rolling_out        TVLQR le long du plan, jusqu'au largage
  phase 2  free_fall          RIEN. u = 0. La balle est balistique et
                              rien de ce que fait le moteur ne change ou
                              elle atterrit.
  phase 3  rolling_in_inside  LQR stationnaire, avec theta replie sur
                              l'angle de garage.

LE RAYON EST DEJA MESURE, ET C'EST CE QUI REND LE PORTAGE POSSIBLE.
control_loop() calcule psi = atan2(u_px - cx, v_px - cy). Les deux memes
composantes donnent hypot(u_px - cx, v_px - cy), donc le rayon en pixels
est disponible sans une ligne de traitement d'image supplementaire. Le
test de phase se fait dessus, exactement comme en simulation.

DETECTION DE PHASE, ET POURQUOI ELLE EST A VERROU
=================================================
En simulation le test est |r - R| < tol ET |r_dot| < tol. Sur le banc
r_dot serait derive d'un rayon bruite a 50 Hz : inutilisable pour une
decision binaire. On utilise donc un VERROU :

  phase 1  tant que t <= Tf
  phase 2  des que t > Tf, jusqu'a ce que r passe sous R_in + tolerance
  phase 3  a partir de la, DEFINITIVEMENT

Une fois la balle rattrapee elle ne repart pas -- et si elle repart,
l'essai est rate de toute facon, donc revenir en phase 2 n'aurait aucune
valeur. Un verrou vaut mieux qu'un test qui oscille sur le bruit.

CE QUI N'EST PAS PORTE, ET IL FAUT LE SAVOIR
============================================
La simulation ne verifie PAS que la balle passe le trou : detect_events.m
et handle_impact.m testent l'azimut du CENTRE de la balle a l'instant ou
elle croise le rayon, pas l'arc balaye par une bille de 25 mm pendant que
le cerceau tourne. Mesure du 11/08 sur le plan courant : il faudrait
128.4 deg d'ouverture pour 90 disponibles, soit -13.8 deg de marge. Un
rattrapage reussi en simulation peut donc etre une COLLISION sur le banc.

Consequence pratique : verifier la marge cote MATLAB
(verification/t4d_feasible_map.m) AVANT de lancer, et commencer par des
essais a vitesse reduite. Ce script ne peut pas le verifier lui-meme --
il ne connait pas la trajectoire du vol, seulement l'etat de largage.
"""

import os

from tasks.t1_config import *          # noqa: F401,F403  (banc, modele, securites)
from tasks import t1_config as _base

# --- Redefinitions propres a T4 ---------------------------------------

# La balle monte a ~120 deg avant de decoller, puis vole. Pendant le vol
# psi n'a plus de sens physique fort (la balle passe pres de l'axe, donc
# l'angle balaie vite), mais il ne doit pas faire sauter la garde.
UNWRAP_PSI     = True
ABORT_PSI_DEG  = 400.0

# CORRIGE LE 28/08 : 25 ETAIT LU SUR LE MAUVAIS CHIFFRE.
# 8 rad/s est la vitesse AU LARGAGE, pas le maximum sur la trajectoire :
# le pompage qui amene la balle a 120 deg monte a 46.30 rad/s sur le plan
# exporte ce jour. A 25, l'anti-windup ecretait la consigne ET recalculait
# u par back-calculation pendant tout le pompage -- le plan aurait ete
# detruit avant meme le largage. 55 laisse 19 % de marge au-dessus du
# plan et reste sous les 251 rad/s auxquels E1 a tourne.
THETA_DOT_MAX  = 55.0

# Le plan exporte le 28/08 monte a 163.91 rad/s^2, donc 140 l'aurait
# ecrete. 200 = 163.91 + 22 % de marge de correction TVLQR, meme
# proportion que les 140 pour 112 en T3. Ce n'est PAS une borne de
# non-glissement -- celle-la est imposee dans le plan, pas ici.
U_SLIP_MAX     = 200.0

# Pose : la balle va moins vite qu'au sommet du looping (7.3 rad/s au
# largage contre 15.6 au pic de T3), mais elle traverse le trou vite.
# 8000 us donne 5.4 mm de file a 7.3 rad/s, soit 0.22 diametre.
EXPOSURE_US    = 8000

# --- Geometrie en pixels, pour le test de phase -----------------------
#
# Le test de phase compare le rayon MESURE, en pixels, aux deux rayons
# d'orbite exportes par MATLAB, en metres. Il faut donc une echelle.
# Elle se deduit du rayon d'orbite exterieur : R_EFF_M metres valent
# HOOP_RADIUS_PX pixels.
#
# HOOP_RADIUS_PX est A MESURER sur une image du banc (rayon du cercle
# decrit par le CENTRE de la balle sur le cerceau exterieur, en pixels
# pleine resolution). Tant qu'il est faux, le test de phase l'est aussi.
R_EFF_M          = 0.09158      # [m] rayon d'orbite exterieur (MATLAB)
HOOP_RADIUS_PX   = 205.2        # [px] MESURE (hoop_geometry_v3, cf. t3_config)

# Pendant le vol la balle peut disparaitre quelques images derriere la
# structure du cerceau interieur. ABORT_LOST_FRAMES = 5 (herite de
# t1_config) avorterait l'essai en plein vol pour 100 ms de perte. 12 =
# 240 ms, encore tres court devant les 155 ms de vol, mais il faut que la
# garde reste utile.
ABORT_LOST_FRAMES = 12

# Estimateur : la DERIVEE FILTREE, pas l'EKF. ekf_v3 porte le modele du
# pendule du cerceau EXTERIEUR (F_N_HZ = 1.3144). Pendant le vol ce modele
# est faux tout court, et en phase 3 la balle est sur la face interne ou
# f_n = 2.2157 Hz. Comme control_loop REMPLACE psi par l'estimation, un
# filtre hors domaine ne degrade pas seulement psi_dot : il corrompt la
# mesure. C'est le meme argument qui interdit --estimator kalman en T3.
DEFAULT_ESTIMATOR = "derivative"

# Bande autour de R_in qui declenche le verrou de phase 3. 4 mm : large
# devant le bruit de detection (0.24-0.63 deg sur l'outer, soit ~1 px)
# et etroit devant les 59 mm qui separent les deux rayons d'orbite.
CATCH_TOL_M      = 4.0e-3

# Fenetre de suivi du detecteur.
#
# 320 ETAIT UNE ERREUR, ET ELLE EST DEJA DOCUMENTEE EN T3 : c'est la
# valeur qui a fait passer la detection de 2.6 ms a 9.8 ms, la boucle
# cessant alors d'etre cadencee par la camera, avec un pas sur cinq a
# 40 ms au lieu de 20. Le chiffre a comparer est le deplacement de la
# balle ENTRE DEUX IMAGES, qui doit tenir dans la DEMI-fenetre :
#   roulage  HOOP_RADIUS_PX * max|psi_dot| / LOOP_HZ
#            = 205.2 * 15.75 / 50 = 64.6 px
#   vol      tangentiel 30.0 px + radial 8.5 px = 38.5 px, moins que le
#            roulage : le vol n'est PAS le cas dimensionnant, contrairement
#            a ce que disait le commentaire precedent.
# 160 donne 80 px de demi-fenetre, soit 24 % de marge sur les 64.6 px.
TRACKING_WIN     = 160

# Fichier de plan exporte par MATLAB (export_t4_for_bench.m).
PLAN_FILE      = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "plans", "T4d_flying_ball_plan.mat")

# Duree de maintien apres le rattrapage, pendant laquelle le LQR de
# phase 3 tient la balle et gare le cerceau.
HOLD_AFTER_CATCH_S = 5.0

# Repertoire propre a T4.
OUTDIR = "data/t4_v3"
