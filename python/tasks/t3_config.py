#!/usr/bin/env python3
"""
t3_config_v3.py -- T3 : looping complet de la balle dans le cerceau exterieur.

Meme principe que t2_config_v3.py : on herite de t1_config_v3 tout ce qui
decrit le BANC (camera, ODrive, cadences, moteur, geometrie) et on ne
redefinit ici que ce qui change parce que la manoeuvre change.

Rien de ce qui touche a l'acquisition n'est modifie : LOOP_HZ, COM_PERIOD,
ENC_PERIOD, ERR_PERIOD, la ROI et le detecteur sont ceux de T1/T2.

CE QUI CHANGE PAR RAPPORT A T1/T2, ET POURQUOI
==============================================

T1 et T2 vivent au voisinage de psi = 0, avec un gain LQR constant. T3 fait
faire un tour complet a la balle en suivant une trajectoire planifiee. Cinq
consequences, dont quatre sont des CORRECTIONS OBLIGATOIRES et non des
reglages :

1. PSI DOIT ETRE DEROULE (unwrap).
   La mesure est psi = atan2(...), donc repliee dans (-pi, pi]. Sur un tour
   complet elle saute de +180 a -180 deg. Ce saut passe ensuite dans le
   derivateur de psi_dot, qui lit un bond de 2*pi en un pas de 20 ms, soit
   ~314 rad/s de vitesse fantome, et dans l'erreur de suivi, qui devient
   fausse de 360 deg. C'est le changement le plus important de T3 et il n'a
   pas d'equivalent dans T1/T2, ou psi ne quitte jamais le voisinage de 0.

2. LA LIMITE D'ABANDON SUR PSI DOIT SAUTER.
   bench_common_v3 coupe des |psi| > ABORT_PSI_DEG. Sur un looping psi va
   jusqu'a -360 deg : la garde se declencherait au premier quart de tour.
   Elle est portee a une valeur qui garde son role (detecter une balle
   perdue ou un decrochage) sans interdire la manoeuvre.

3. THETA_DOT_MAX DOIT ETRE RELEVE.
   La trajectoire T3c/F demande 25.04 rad/s au cerceau, contre les 20 rad/s
   de T1/T2. L'anti-windup de la boucle ecreterait la consigne ET
   recalculerait u par back-calculation : le plan serait detruit exactement
   la ou il compte. Le moteur n'est pas en cause -- l'identification E1 a
   tourne avec vel_limit a 40 tr/s = 251 rad/s, et a 100 rad/s le couple
   visqueux ne represente que 1.1 % de TAU_MAX.

4. U_SLIP_MAX DOIT ETRE REVU.
   Le plan T3c/F monte a max|u| = 112.24 rad/s^2, au-dessus des 100 de
   T1/T2 : la saturation ecreterait le plan. Surtout, la constante 100 n'est
   valable qu'au fond du cerceau. La vraie borne de non-glissement est
   mu*N-dependante et s'effondre au sommet, ou N tend vers zero -- c'est
   precisement la que le plan passe. Le planificateur MATLAB impose deja
   cette borne exacte a chaque noeud (enforce_no_slip, dynamics/
   contact_forces.m), donc la trajectoire EST admissible ; la constante
   ici ne doit plus servir qu'a rattraper une divergence, pas a filtrer le
   plan. Elle est relevee juste au-dessus du plan.

   COROLLAIRE A ASSUMER : a 140, la saturation ne protege PLUS du
   glissement. Le plan monte a 112, donc le TVLQR peut ajouter jusqu'a
   28 rad/s^2 de correction -- et il le fera au sommet, la ou N tend vers
   zero et ou la vraie borne mu*N est la plus basse. Ce n'est pas
   dangereux pour le banc (le moteur dispose de 483 rad/s^2), mais un
   essai ou la balle glisse est HORS MODELE et donc ininterpretable. Le
   critere a posteriori est un decrochage soudain de psi par rapport au
   plan au voisinage du sommet : il se lit dans le journal, pas en ligne.

5. LE GAIN EST VARIANT DANS LE TEMPS.
   T1/T2 utilisent un K constant. T3 utilise K(t) issu du TVLQR. Cela
   demande la seule modification de bench_common_v3.py que T3 impose --
   voir l'en-tete de t3_run_v3.py.

RISQUE MATERIEL A VERIFIER AVANT DE LANCER
==========================================
Au sommet du looping la balle doit passer a psi_dot >= sqrt(g/R_eff) =
10.35 rad/s, soit 0.95 m/s au centre de la balle. Avec EXPOSURE_US = 18000
elle se deplace de 17 mm pendant la pose, pour un diametre de 25 mm : le
file est de 0.7 diametre et la detection peut echouer LA OU C'EST LE PLUS
CRITIQUE. EXPOSURE_US est donc reduit ici, ce qui coute de la lumiere : il
faut verifier la detection avec ball_detection.py --preview avant l'essai,
et remonter GAIN_CAM ou l'eclairage si le blob devient trop sombre.
"""

import math
import os

from tasks.t1_config import *          # noqa: F401,F403  (banc, modele, securites)
from tasks import t1_config as _base

# --- Herite tel quel : banc, camera, ODrive, cadences, moteur, modele ---
#
# L'import est en * comme dans t2_config_v3, et PAS en liste explicite.
# La version precedente listait a la main une trentaine de constantes et en
# oubliait dix-huit : THRESHOLD, DOWNSAMPLE, TRACKING_WIN, BALL_SIZE,
# ROI_MARGIN, HOOP_CENTRE_PX, USE_HOUGH_CENTRE, CAM_LATENCY_S,
# SIGMA_PSI_DEG, Q_ACCEL, RAMP_DOWN_S, GAINS et les six ARM_*.
# Consequence : AttributeError sur cfg.THRESHOLD des bc.setup_camera(), donc
# plantage au premier essai reel -- et --dry-run ne le voyait pas, puisqu'il
# sort avant de toucher au materiel. Une liste explicite est une source
# d'erreur silencieuse : ce qui est redefini pour T3 est de toute facon
# visible plus bas, c'est la seule chose qui compte.


# --- Redefinitions propres a T3 (voir l'en-tete pour le pourquoi) -------

# (1)+(2) La balle fait un tour complet : psi deroule descend a -360 deg.
# La garde garde son role -- detecter une divergence ou une detection
# aberrante -- mais a un seuil que la manoeuvre ne franchit pas. 600 deg
# laisse deux tiers de revolution de depassement au-dela de la cible, la
# simulation en boucle fermee sur T3b ayant mesure un overshoot a -443 deg.
UNWRAP_PSI     = True
ABORT_PSI_DEG  = 600.0

# (3) La trajectoire demande 25.04 rad/s ; 35 laisse de la marge de suivi
# sans approcher quoi que ce soit de physique. vel_limit ODrive suivra a
# 1.2 x cette valeur = 42 rad/s = 6.7 tr/s, tres en dessous des 40 tr/s
# utilises en E1.
THETA_DOT_MAX  = 35.0

# (4) LE PLAN COULOMB DEMANDE 247 rad/s^2, PAS 112.
# Valeur d'origine : 140, dimensionnee sur le plan T3c/F qui montait a
# 112.24 rad/s^2. Le plan T3_loop_the_loop_coulomb (export du 28/08 :
# Tf = 1.4134 s, 60 noeuds) monte a 247.03. A 140, t3_run_v3.py refuse de
# partir (check_plan_against_config) ; s'il partait, la saturation
# ecreterait le plan sur toute la partie qui compte et le test ne voudrait
# rien dire.
#
# 300 = 247 + 21 % de marge de correction TVLQR, la meme proportion que
# 140 pour 112. LE COUPLE SUIT : a 18.85 rad/s la borne moteur vaut
# (TAU_MAX - B_MOTEUR * theta_dot) / I_TOTAL = 483.6 rad/s^2, donc 300 ne
# demande que 62 % de ce dont le banc dispose.
#
# ATTENTION AU DESACCORD DE TAU_MAX ENTRE LES DEUX COTES. Le plan est
# resolu avec u_max = cfg.tau_max / total_inertia cote MATLAB, et
# sim_config.m y porte tau_max = 0.310 N.m (soit 10 A a Kt = 0.031) alors
# que ce fichier porte TAU_MAX = 0.5425 (17.5 A). Le planificateur est
# donc le plus conservateur des deux et le plan tient dans le budget du
# banc quel que soit celui qui a raison -- mais les deux chiffres ne
# peuvent pas etre vrais en meme temps, et le memoire les cite tous les
# deux. A trancher en lisant la limite reellement configuree dans l'ODrive.
#
# COROLLAIRE INCHANGE ET AGGRAVE : a 300 la saturation ne protege plus du
# tout du glissement, et la marge de non-glissement du plan Coulomb est
# serree PAR CONSTRUCTION -- c'est ce que la contrainte impose, elle
# sature. Le critere a posteriori reste le meme : decrochage soudain de
# psi par rapport au plan au voisinage du sommet, lu dans le journal.
U_SLIP_MAX     = 300.0

# (5) Pose reduite : voir "RISQUE MATERIEL" en tete de fichier.
EXPOSURE_US    = 9000
# Seuil valide avec l'eclairage et les reglages de gain du banc.
THRESHOLD      = 90
# ANALOGUE_GAIN etait du code mort : open_camera lit GAIN_CAM, herite de
# t1_config_v3 par l'import *. Regler le gain camera pour T3 se fait en
# redefinissant GAIN_CAM ici, pas ANALOGUE_GAIN.

# (5b) LA POSE N'EST PAS LE SEUL PROBLEME DE DETECTION, NI LE PRINCIPAL.
# Le file pendant la pose vaut 5.7 mm a 6000 us, soit 0.23 diametre : c'est
# regle. Mais ENTRE DEUX IMAGES la balle parcourt R_eff * psi_dot * dt =
# 0.0916 * 10.35 * 0.02 = 19 mm, soit 0.76 diametre. En pixels cela fait
# HOOP_RADIUS_PX * psi_dot * dt, a comparer a TRACKING_WIN / 2 : si le
# deplacement sort de la fenetre de suivi, le detecteur perd la balle au
# sommet, exactement la ou la manoeuvre se joue.
#
# HOOP_RADIUS_PX est A MESURER sur une image du banc (rayon du cercle
# decrit par le CENTRE de la balle, en pixels pleine resolution). Tant
# qu'il vaut None, t3_run_v3.py refuse de partir : c'est une verification
# geometrique, pas un reglage a deviner.
R_EFF_M          = 0.09158      # [m] rayon effectif de la trajectoire
HOOP_RADIUS_PX   = 208.0        # [px] A MESURER, puis figer ici

# Le deplacement inter-image au pic vaut HOOP_RADIUS_PX * psi_dot * dt =
# 208 * 15.55 / 50 = 65 px : la demi-fenetre doit le depasser, pas plus.
# 320 le depassait d'un facteur 2.5 et coutait 25 fois la surface de la
# fenetre T1/T2 (64), d'ou une detection a 9.8 ms contre 2.6 ms en T2 sur
# une ROI POURTANT PLUS GRANDE (565x559 contre 525x525). La boucle cessait
# alors d'etre cadencee par la camera -- 'attente image' tombe de 55 % en T2
# a 20 % en T3 -- et un pas sur cinq durait 40 ms au lieu de 20. C'est cette
# distribution bimodale de dt, et non un decalage moyen, qui explique
# l'apparente contradiction entre 'TOTAL 21 ms -> 47.6 Hz' et 'cadence
# mediane 52 Hz' dans les resumes de fin d'essai.
#
# Consequence directe sur l'estimation : ekf_v3 compense le retard camera
# sur n_delay = round(CAM_LATENCY_S / DT_NOMINAL) = 1 pas FIXE. Sur un pas
# de 40 ms la compensation porte sur le mauvais horizon, et 20 ms d'erreur
# a 15 rad/s font 17 deg d'erreur sur psi -- du meme ordre que ce que le
# TVLQR cherche a corriger. En T1/T2 (psi_dot < 1 rad/s) c'etait invisible.
#
# 160 laisse 80 px de demi-fenetre, soit 23 % de marge sur les 65 px du pic.
#
# RECALCUL POUR LE PLAN COULOMB (28/08). Il fait le meme tour en 1.4134 s
# au lieu de 2.38, donc plus vite : max |psi_dot| = 16.68 rad/s au lieu de
# 15.55. Le deplacement inter-image devient HOOP_RADIUS_PX * psi_dot * dt
# = 205.2 * 16.68 / 50 = 68.5 px, contre 80 px de demi-fenetre : il reste
# 14 % de marge au lieu de 23. Ca passe, et c'est la seule raison pour
# laquelle TRACKING_WIN n'est pas releve -- le cout de la fenetre est
# MESURE (9.8 ms de detection a 320, contre 2.6 ms en T2), son benefice
# ici ne le serait pas. Si les journaux du premier essai Coulomb montrent
# des images perdues au voisinage du sommet, c'est ce nombre qu'il faut
# monter, a 192 (96 px, 40 % de marge) avant tout autre reglage.
#
# Le file pendant la pose suit la meme proportion : R_eff * psi_dot * t_exp
# = 0.09158 * 16.68 * 9e-3 = 13.7 mm, soit 0.55 diametre a EXPOSURE_US =
# 9000. C'etait 12.8 mm avec le plan precedent.
TRACKING_WIN     = 160

# Rayon qui DIMENSIONNE LA ROI, distinct de HOOP_RADIUS_PX ci-dessus qui
# decrit l'orbite du CENTRE. Il faut au minimum 208 + 28 = 236 px pour que
# la balle entiere tienne dans la ROI ; 250 avec ROI_MARGIN = 0.05 donne
# une demi-ROI de 262 px, soit 26 px de garde au-dela du bord de la balle.
# Voir roi_from_configured_hoop() pour la mesure du 10/08 qui l'impose.
#
# ATTENTION, CECI NE CORRIGE QUE LE ROGNAGE, PAS LE BRUIT. La ROI grandit
# (437 -> 525 px de cote), et le bruit sur psi croit avec elle : mesure au
# 10/08, l'outer a 565x559 donne un ecart-type de 0.24-0.63 deg contre
# 0.05 sur l'inner a 243x243, un facteur 10 du seul fait de la taille. Les
# deux exigences tirent en sens oppose et seule une ROI EN COURONNE
# (orbite +/- ~40 px, sans le grand vide central) les satisfait ensemble.
# Chantier separe : il change le chemin de donnees du detecteur et
# toucherait aussi la chaine de mesure de T1.
#ROI_RADIUS_PX    = 250.0

HOOP_CENTRE_PX = (377.92, 271.60)
HOOP_RADIUS_PX = 205.2
ROI_RADIUS_PX = 230

# Retard de la boucle de vitesse, mesure par vel_loop_id_v3.py.
# A RENSEIGNER APRES MESURE. Tant qu'il vaut 0.0, control_loop se comporte
# exactement comme en T1/T2 : le retour d'etat regule sur la consigne
# integree et non sur le cerceau reel, ce que les journaux T3 invalident
# (|vel_cmd - vel_enc| a 1.5-2.7 rad/s RMS contre 0.15-0.40 en T1/T2).
TAU_VEL_S      = 0.0

# (5c) CAM_LATENCY_S = 19.3 ms a ete mesuree A 18000 us DE POSE. La latence
# contient environ la moitie du temps d'integration, donc a 9000 us elle
# devrait descendre d'a peu pres 5 ms -- raisonnement, pas mesure, d'ou
# l'absence de valeur de remplacement ici.
#
# CETTE CONSTANTE EST ACTIVE depuis que --estimator ekf est le defaut
# (t3_run_v3.py) : ekf_v3.build_estimator compense sur
#     n_delay = int(round(CAM_LATENCY_S / DT_NOMINAL)).
# Elle n'etait que journalisee tant que T3 tournait sur "derivative", qui
# renvoie n_delay = 0.
#
# L'INEXACTITUDE N'EST PAS BLOQUANTE, et c'est le calcul qui le dit plutot
# qu'un avis : n_delay est arrondi au pas entier, donc a DT_NOMINAL = 20 ms
# TOUTE latence comprise entre environ 11 et 29 ms donne n_delay = 1. Que
# la vraie valeur soit 19.3 ou 14 ms ne change donc rien a la compensation
# appliquee. La re-mesure (latency_bench.py, A 9000 us) reste souhaitable
# pour pouvoir CITER un chiffre dans le memoire ; elle ne conditionne pas
# l'usage de l'EKF.
CAM_LATENCY_S  = _base.CAM_LATENCY_S

# (5d) L'ecart entre la consigne de vitesse et la vitesse reelle du cerceau
# est LA quantite qui dit si l'hypothese de boucle interne parfaite tient a
# 25 rad/s et 112 rad/s^2. A ENC_PERIOD = 10 elle n'est vue qu'a 5 Hz, soit
# une dizaine de points sur tout le looping. 2 la porte a 25 Hz pour deux
# aller-retours USB par tour, ce que le budget supporte.
#
# (5e) LE RETOUR D'ETAT UTILISE MAINTENANT L'ENCODEUR.
# =====================================================
# MESURE sur les huit essais du 10/08 (Analyse/t3_delay_budget.py) : l'ecart
# RMS entre le theta REGULE (consigne integree) et le theta ENCODEUR vaut
#     t01 11.1   t02  5.4   t03  9.8   t06 8.3   t07 11.7   t08 7.7  deg
#     t04 41.6   t05 238.9 deg
# et les deux valeurs hautes appartiennent a des essais qui divergent. Avec
# USE_MEASURED_THETA = False cet ecart est une erreur que le regulateur NE
# VOIT PAS et qui entre pourtant dans u par K1 et K2.
#
# Le sens de causalite n'est pas etabli : la saturation peut creuser l'ecart
# autant que l'ecart peut causer la divergence. Mais 238.9 deg ne s'explique
# pas par les 5 % de sous-livraison du variateur (mesure : vel_enc =
# 0.95 * vel_cmd, retard d'actionnement 0.0-0.9 ms), et le correctif se teste
# en un essai.
#
# ENC_PERIOD = 1 N'EST PAS UN CHOIX : validate_config() refuse de demarrer
# avec USE_MEASURED_THETA = True et ENC_PERIOD != 1, puisque l'encodeur
# entre alors dans la loi de commande et doit etre lu a chaque tour.
#
# RISQUE A VERIFIER DES LE PREMIER ESSAI. Passer de 5 a 1 quintuple le
# trafic USB de la voie encodeur, et la note (5d) ci-dessus ne donne le
# budget comme sur qu'a partir de 2. C'est le mecanisme qui avait bride la
# boucle T1 a ~36 Hz (lecture d'Iq a chaque iteration) et fait diverger cinq
# essais sur cinq. TEST D'ACCEPTATION : sur le premier essai, verifier que
# le dt median reste a ~20 ms et que la fraction de pas > 1.5 x median ne
# depasse pas les 2.8 % mesures au 10/08. Si la cadence s'effondre, revenir
# a USE_MEASURED_THETA = False / ENC_PERIOD = 5 (sauvegarde dans
# backup_avant_modifs_T3/) : une boucle lente est un mal pire que le biais
# de theta qu'on cherche a corriger.
USE_MEASURED_THETA = True
ENC_PERIOD     = 1

# Dans t3_config_v3.py
#ENC_PERIOD = 5   # 10 Hz -- valeur d'avant (5e), incompatible avec
#                 # USE_MEASURED_THETA = True
#ERR_PERIOD = 10  # 5 Hz
#IQ_PERIOD  = 0   # déjà désactivé

# Fichier de plan exporte par MATLAB (export_tvlqr_for_bench.m).
#
# DEFAUT CHANGE LE 28/08 : le plan Coulomb remplace T3c/F.
#   T3c_loop_no_slip_plan.mat     Tf = 2.38 s, max |u| = 112.24 rad/s^2
#   T3_loop_the_loop_coulomb_...  Tf = 1.4134 s, max |u| = 247.03 rad/s^2
# Le second impose la contrainte de Coulomb noeud par noeud
# (enforce_no_slip, contact_margin = 0.20 g) et fait le tour presque deux
# fois plus vite. Il exige U_SLIP_MAX = 300, voir (4) ci-dessus.
#
# L'ancien plan reste dans plans/ et se rejoue sans rien modifier ici :
#     python3 t3_run_v3.py --plan plans/T3c_loop_no_slip_plan.mat
# mais U_SLIP_MAX est alors trois fois plus haut que ce que ce plan-la
# demande, donc la saturation ne le borne plus : c'est le seul point ou
# les deux plans ne sont pas interchangeables a config egale.
PLAN_FILE      = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "plans", "T3_loop_the_loop_coulomb_plan.mat")

# Duree de maintien apres la fin de la trajectoire. Au-dela de Tf le TVLQR
# tient son dernier gain comme regulateur stationnaire, exactement comme
# tvlqr_controller.m cote simulation.
HOLD_AFTER_TF_S = 4.0

# Repertoire propre a T3 : _base.OUTDIR vaut "data/t1_v3" et melanger les
# campagnes rendrait l'analyse ambigue.
OUTDIR = "data/t3_v3"
