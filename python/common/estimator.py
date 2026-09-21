"""
estimator_v3.py -- Estimation de (psi, psi_dot) avec compensation de retard.

Reprend l'approche de Gurtner & Zemanek (IFAC 2017, section 4.2), qui
differe sur trois points de la premiere version de t1_run.py.

1. Filtre de Kalman au lieu d'une derivee filtree
-------------------------------------------------
Une dérivée suivie d'un passe bas ne connait rien de la physique du système.
Elle doit choisir bruit vs retard de phase. Le filtre de Kalman utilise le 
modèle linéaire de la balle, donc il sait à quoi ressemble une trajectoire 
plausible. 

Il extrait ensuite psi_dot avec beaucoup moins de retard pour un bruit égal.
AA4CC utilise un EKF. 

Autour de psi = 0, le modèle est déjà linéaire à quelques % près. Un filtre
linéaire suffit donc, et évite d'utiliser une jacobienne pour les calculs. 

2. Compensation du retard de mesure par prédiction
--------------------------------------------------
L'image rendue par la caméra a été exposée CAM_LATENCY_S plus tôt (~19 ms).
La mesure décrit donc l'état de la balle à un instant antérieur, pas actuel.

Ca introduit un retard pur dans la boucle, ce qui réduit la marge de phase, 
très fortement pour un système aussi peu amorti, qui plus est. La parade 
(aussi utilisée par AA4CC) est de maintenir l'estimation à l'instant t-1, puis
la propager à l'avant de L pas avec l'historique des commandes déjà appliquées, 
qui lui est connu exactement.

    x_hat[k-L]  <- correction par la mesure psi[k]
    x_hat[k-L+1] = f(x_hat[k-L],   u[k-L])
    ...
    x_hat[k]     = f(x_hat[k-1],   u[k-1])

3. theta et theta_dot par integration de la commande
-----------------------------------------------------
AA4CC n'utilise pas l'encodeur pour theta, mais l'entrée u intégrée. C'est 
cohérent avec l'hypothèse de boucle interne parfaite sur laquelle le LQR est
créé. Lire l'encodeur donne les valeurs réelles, qui trainent derrière si la 
boucle interne est lente. Le retour d'état réagit à cet écart là et referme 
une boucle qui parasite complètement la dynamique. 
Pour trancher sur une mesure plutôt que le principe, les mesures sont collectées
à titre de comparaison.

Modele utilisé
--------------
Autour de psi = 0, avec u = theta_ddot :

    d/dt [psi    ]   [  0     1  ] [psi    ]   [ 0  ]
         [psi_dot] = [ a21   a22 ] [psi_dot] + [ b2 ] u

    a21 = -omega_n^2          = -68.20   (f_n = 1.3144 Hz, mesuree en E2)
    a22 = -2 zeta omega_n     = -0.545   (zeta = 0.033, decrement E2)
    b2                        = -0.402

Les coeffs ne sont pas redéfinis mais donnés dans la configuration. Le signe de 
b2 est négatif, car la convention MATLAB et la convention Python sont inversées 
l'une par rapport à l'autre! 

[Gurtner, M. & Zemanek, J., "Ball in double hoop: demonstration
model for numerical optimal control", IFAC-PapersOnLine 50(1), 2017]
"""

import math
import numpy as np


class DelayCompensatedKF:
    """Kalman lineaire sur (psi, psi_dot) + prediction de L pas."""

    def __init__(self, a21, a22, b2, sigma_psi, q_accel, n_delay):
        self.a21, self.a22, self.b2 = a21, a22, b2
        self.R = np.array([[sigma_psi ** 2]])       # variance de mesure
        self.q_accel = q_accel                      # densite de bruit modele
        self.n_delay = int(n_delay)
        self.H = np.array([[1.0, 0.0]])

        self.x = np.zeros(2)                        # etat a l'instant retarde
        self.P = np.diag([sigma_psi ** 2, 1.0])

        # Historique des commandes. INVARIANT (L = n_delay) : a l'entree de
        # update() qui traite la mesure psi[k],
        #     u_hist = [u[k-L-1], u[k-L], ..., u[k-1]]      (L+1 elements)
        # donc u_hist[0] est la commande a appliquer pour propager l'etat
        # retarde d'un pas AVANT la correction, et u_hist[1:] sont les L
        # commandes deja envoyees qui servent a la prediction en avant.
        # La v3 initiale dimensionnait cette liste a max(L, 1) : le meme
        # echantillon servait alors aux deux roles et la propagation etait
        # decalee d'un pas.
        self.u_hist = [0.0] * (self.n_delay + 1)
        self.ready = False
        self._fg_cache = {}

    # -- Modele discret ----------------------------------------------------
    def _FG(self, dt):
        """Discretisation bilineaire (Tustin).

        Euler explicite donne pour cet oscillateur des valeurs propres de
        module |1 +/- j w_n dt| = 1.014 a dt = 20 ms : la discretisation
        elle-meme est legerement instable, ce que la correction masque sans
        le supprimer. Tustin conserve le demi-plan gauche pour tout dt et
        coute une inversion 2x2, negligeable et mise en cache puisque dt
        varie peu.
        """
        key = round(dt, 9)
        cached = self._fg_cache.get(key)
        if cached is not None:
            return cached
        A = np.array([[0.0, 1.0], [self.a21, self.a22]])
        B = np.array([0.0, self.b2])
        M = np.linalg.inv(np.eye(2) - 0.5 * dt * A)
        F = M @ (np.eye(2) + 0.5 * dt * A)
        G = M @ (B * dt)
        if len(self._fg_cache) > 256:
            self._fg_cache.clear()
        self._fg_cache[key] = (F, G)
        return F, G

    def _Q(self, dt):
        """Bruit de modele vu comme une acceleration aleatoire : rend compte
        du frottement mal modelise et des perturbations de contact."""
        return self.q_accel * np.array([[dt ** 3 / 3.0, dt ** 2 / 2.0],
                                        [dt ** 2 / 2.0, dt]])

    def _step(self, x, P, u, dt):
        F, G = self._FG(dt)
        return F @ x + G * u, F @ P @ F.T + self._Q(dt)

    # -- Interface publique ------------------------------------------------
    def reset(self, psi0):
        self.x = np.array([psi0, 0.0])
        self.P = np.diag([self.R[0, 0], 1.0])
        self.u_hist = [0.0] * (self.n_delay + 1)
        self.ready = True
        self.x_now = self.x.copy()

    def update(self, psi_meas, dt):
        """Corrige l'etat retarde avec la mesure, puis le propage jusqu'a
        maintenant. Retourne l'estimation COURANTE (psi, psi_dot)."""
        if not self.ready:
            self.reset(psi_meas)
            return self.x.copy()

        # --- Propagation d'un pas a l'instant retarde, puis correction ---
        # u_hist[0] = u[k-L-1] : la commande qui a agi entre l'instant
        # retarde precedent et celui que la mesure psi[k] decrit.
        self.x, self.P = self._step(self.x, self.P, self.u_hist[0], dt)

        S = self.H @ self.P @ self.H.T + self.R
        Kg = (self.P @ self.H.T) @ np.linalg.inv(S)
        # float() sur un tableau de forme (1,) est une erreur depuis
        # NumPy 2.0 (et une DeprecationWarning avant) : on indexe.
        innov = float(psi_meas) - float((self.H @ self.x)[0])
        self.x = self.x + (Kg.flatten() * innov)
        self.P = (np.eye(2) - Kg @ self.H) @ self.P

        # --- Prediction en avant sur les commandes deja appliquees ---
        # u_hist[1:] = u[k-L] ... u[k-1], soit exactement les L commandes
        # envoyees depuis l'instant que la mesure decrit.
        x_now = self.x.copy()
        P_now = self.P.copy()
        for u_k in self.u_hist[1:]:
            x_now, P_now = self._step(x_now, P_now, u_k, dt)

        self.x_now = x_now
        return x_now

    def push_command(self, u):
        """A appeler apres avoir envoye u au moteur, pour que la prediction
        du tour suivant dispose de l'historique correct.

        DOIT etre appelee a CHAQUE pas de boucle, y compris sur les images
        ou la balle n'est pas detectee (pousser 0.0, la commande y etant
        gelee). Un appel manquant desynchronise l'historique et la
        compensation de retard porte alors sur les mauvaises commandes.
        """
        self.u_hist.append(float(u))
        while len(self.u_hist) > self.n_delay + 1:
            self.u_hist.pop(0)

    def predict_open_loop(self, dt):
        """Propagation sans mesure, quand la balle n'est pas detectee.
        Preferable a un gel de l'estimation sur une ou deux images.

        Meme decoupage que update() : un pas sur l'etat retarde avec
        u[k-L-1], puis L pas de prediction. La covariance croit, ce qui est
        le comportement voulu : l'estimation se degrade tant qu'aucune
        mesure ne la recale.
        """
        if not self.ready:
            return self.x.copy()
        self.x, self.P = self._step(self.x, self.P, self.u_hist[0], dt)
        x_now, P_now = self.x.copy(), self.P.copy()
        for u_k in self.u_hist[1:]:
            x_now, P_now = self._step(x_now, P_now, u_k, dt)
        self.x_now = x_now
        return x_now


def build_estimator(cfg, dt_nominal):
    """Construit le filtre a partir de t1_config.

    DOMAINE DE VALIDITE. Le modele est linearise autour de psi = 0 : le
    rappel vaut A21 * psi au lieu de A21 * sin(psi). L'ecart est de 1 % a
    14 deg, de 10 % a 45 deg, et le terme change carrement de signe au-dela
    de 180 deg. Le filtre n'est donc PAS utilisable sur une manoeuvre qui
    quitte le voisinage du fond -- typiquement T3, ou psi descend a
    -360 deg. Le refus est ici plutot que dans le script d'essai, pour que
    tout appelant en beneficie.
    """
    if getattr(cfg, "UNWRAP_PSI", False):
        raise SystemExit(
            "estimator_v3 : le filtre de Kalman est lineaire autour de "
            "psi = 0 ;\nla configuration a UNWRAP_PSI = True, donc psi "
            "parcourt un tour complet et\nle modele y est faux (le rappel "
            "change de signe des |psi| > 90 deg).\nDe plus control_loop() "
            "remplace psi par l'estimation, donc un filtre hors\ndomaine "
            "corrompt la MESURE et pas seulement psi_dot.\n"
            "Utiliser --estimator derivative pour cette campagne.")

    n_delay = int(round(cfg.CAM_LATENCY_S / dt_nominal))
    kf = DelayCompensatedKF(
        a21=cfg.A21, a22=cfg.A22, b2=cfg.B2,
        sigma_psi=math.radians(cfg.SIGMA_PSI_DEG),
        q_accel=cfg.Q_ACCEL,
        n_delay=n_delay)
    return kf, n_delay
