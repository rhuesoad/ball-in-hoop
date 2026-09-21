"""
ekf_v3.py -- Estimation de (psi, psi_dot) par filtre de Kalman ETENDU.

Meme role que estimator_v3.DelayCompensatedKF, meme interface publique, meme
compensation de retard. La seule difference est le modele porte par le filtre,
et c'est cette difference qui autorise T3.

POURQUOI UN EKF ET PAS LE FILTRE LINEAIRE
=========================================
estimator_v3 porte le modele linearise autour de psi = 0 :

    psi_ddot = A21 psi + A22 psi_dot + B2 u

alors que la dynamique reelle a un rappel gravitaire en sinus :

    psi_ddot = A21 sin(psi) + A22 psi_dot + B2 u

L'ecart vaut 1 % a 14 deg, 10 % a 45 deg, et A21*psi change de signe par
rapport a A21*sin(psi) des |psi| > 180 deg -- c'est-a-dire AU SOMMET du
looping, la ou la manoeuvre se joue. Aucun reglage de Q ou de R ne rattrape
un modele dont le signe est faux, d'ou le refus porte par
estimator_v3.build_estimator des que UNWRAP_PSI est vrai.

L'EKF garde la structure du filtre (prediction, innovation, gain de Kalman)
mais relinearise la jacobienne a chaque pas autour de l'estimation courante :

    F_c(psi) = [    0            1   ]
               [ A21 cos(psi)   A22  ]

En psi = 0 on retrouve exactement estimator_v3 : cos(0) = 1. Le filtre
lineaire est donc le cas particulier de celui-ci, ce qui rend les deux
comparables sur T1 et T2 sans changer autre chose que le nom de la classe.

C'est l'approche de AA4CC [Gurtner & Zemanek, IFAC-PapersOnLine 50(1), 2017,
sect. 4.2], qui utilise un EKF et non un filtre lineaire.

CE QUI RESTE VALABLE DE estimator_v3
====================================
1. Compensation du retard de mesure par prediction. L'image a ete exposee
   CAM_LATENCY_S plus tot : l'etat est corrige a l'instant retarde puis
   propage jusqu'a maintenant avec l'historique EXACT des commandes deja
   appliquees. Inchange.
2. Propagation en boucle ouverte quand la balle n'est pas detectee, avec
   croissance de la covariance. Inchange.
3. Discretisation bilineaire de la covariance. Euler explicite donne pour cet
   oscillateur des valeurs propres de module > 1 a dt = 20 ms ; Tustin
   conserve le demi-plan gauche pour tout dt. Ici la jacobienne depend de
   l'etat, donc le cache de estimator_v3 (indexe sur dt seul) ne tient plus
   et l'inversion 2x2 est refaite a chaque pas -- elle coute quelques
   microsecondes, sans commune mesure avec les 13 ms d'attente image.

CE QUI CHANGE
=============
L'etat est propage par RK4 sur la dynamique NON LINEAIRE, pas par F @ x.
Sur un pas de 20 ms et une pulsation propre de 8.3 rad/s (outer) ou
13.9 rad/s (inner), Euler explicite sur la trajectoire elle-meme accumule
une erreur visible sur les L pas de prediction ; RK4 la rend negligeable
pour un cout de quatre evaluations d'un modele a deux etats.

psi n'est PAS ramene dans (-pi, pi]. Le filtre recoit le psi DEROULE de
control_loop (UNWRAP_PSI) et doit le rendre deroule, sinon l'erreur de suivi
saute de 360 deg. sin() et cos() s'en accommodent sans precaution ; c'est
l'innovation qui l'exigerait si la mesure etait repliee, et elle ne l'est pas.
"""

import math
import numpy as np


class ExtendedKF:
    """Kalman etendu sur (psi, psi_dot) + prediction de L pas.

    Interface identique a estimator_v3.DelayCompensatedKF : reset, update,
    push_command, predict_open_loop. Les deux sont donc interchangeables
    dans control_loop sans y toucher.
    """

    def __init__(self, a21, a22, b2, sigma_psi, q_accel, n_delay):
        self.a21, self.a22, self.b2 = a21, a22, b2
        self.R = np.array([[sigma_psi ** 2]])       # variance de mesure
        self.q_accel = q_accel                      # densite de bruit modele
        self.n_delay = int(n_delay)
        self.H = np.array([[1.0, 0.0]])

        self.x = np.zeros(2)                        # etat a l'instant retarde
        self.P = np.diag([sigma_psi ** 2, 1.0])
        self.x_now = self.x.copy()

        # Historique des commandes. INVARIANT (L = n_delay) : a l'entree de
        # update() qui traite la mesure psi[k],
        #     u_hist = [u[k-L-1], u[k-L], ..., u[k-1]]      (L+1 elements)
        # u_hist[0] propage l'etat retarde d'un pas AVANT la correction,
        # u_hist[1:] sont les L commandes deja envoyees qui servent a la
        # prediction en avant.
        self.u_hist = [0.0] * (self.n_delay + 1)
        self.ready = False

    # -- Modele -----------------------------------------------------------
    def _f(self, x, u):
        """Dynamique non lineaire. C'est ICI et nulle part ailleurs que
        l'EKF differe du filtre lineaire : sin(psi) au lieu de psi."""
        return np.array([x[1],
                         self.a21 * math.sin(x[0]) + self.a22 * x[1]
                         + self.b2 * u])

    def _jacobian(self, x):
        """Jacobienne de _f par rapport a l'etat, evaluee en x.

        d/dpsi [A21 sin(psi)] = A21 cos(psi). En psi = 0 cela vaut A21 et on
        retombe sur estimator_v3 ; a 90 deg le rappel disparait ; au-dela de
        180 deg il change de signe, ce qui est precisement le comportement
        que le filtre lineaire ne pouvait pas representer.
        """
        return np.array([[0.0, 1.0],
                         [self.a21 * math.cos(x[0]), self.a22]])

    def _Q(self, dt):
        """Bruit de modele vu comme une acceleration aleatoire : rend compte
        du frottement mal modelise et des perturbations de contact."""
        return self.q_accel * np.array([[dt ** 3 / 3.0, dt ** 2 / 2.0],
                                        [dt ** 2 / 2.0, dt]])

    def _step(self, x, P, u, dt):
        """Un pas : etat par RK4 sur la dynamique non lineaire, covariance
        par la jacobienne discretisee en bilineaire."""
        k1 = self._f(x, u)
        k2 = self._f(x + 0.5 * dt * k1, u)
        k3 = self._f(x + 0.5 * dt * k2, u)
        k4 = self._f(x + dt * k3, u)
        x_new = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        # Jacobienne au milieu du pas : centrer l'evaluation donne l'ordre 2
        # en dt sur la propagation de covariance, pour le meme cout.
        Fc = self._jacobian(0.5 * (x + x_new))
        M = np.linalg.inv(np.eye(2) - 0.5 * dt * Fc)
        Fd = M @ (np.eye(2) + 0.5 * dt * Fc)
        P_new = Fd @ P @ Fd.T + self._Q(dt)
        return x_new, P_new

    # -- Interface publique -----------------------------------------------
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
        self.x, self.P = self._step(self.x, self.P, self.u_hist[0], dt)

        S = self.H @ self.P @ self.H.T + self.R
        Kg = (self.P @ self.H.T) @ np.linalg.inv(S)
        innov = float(psi_meas) - float((self.H @ self.x)[0])
        self.x = self.x + (Kg.flatten() * innov)
        # Forme de Joseph : (I-KH)P(I-KH)' + KRK'. Plus couteuse que
        # (I-KH)P, mais elle garde P symetrique definie positive meme quand
        # la jacobienne varie fortement d'un pas a l'autre -- ce qui est le
        # cas au passage du sommet, ou A21 cos(psi) traverse zero.
        IKH = np.eye(2) - Kg @ self.H
        self.P = IKH @ self.P @ IKH.T + Kg @ self.R @ Kg.T

        # --- Prediction en avant sur les commandes deja appliquees ---
        x_now, P_now = self.x.copy(), self.P.copy()
        for u_k in self.u_hist[1:]:
            x_now, P_now = self._step(x_now, P_now, u_k, dt)

        self.x_now = x_now
        return x_now

    def push_command(self, u):
        """A appeler apres avoir envoye u au moteur, pour que la prediction
        du tour suivant dispose de l'historique correct.

        DOIT etre appelee a CHAQUE pas de boucle, y compris sur les images ou
        la balle n'est pas detectee (pousser 0.0, la commande y etant gelee).
        Un appel manquant desynchronise l'historique et la compensation de
        retard porte alors sur les mauvaises commandes.
        """
        self.u_hist.append(float(u))
        while len(self.u_hist) > self.n_delay + 1:
            self.u_hist.pop(0)

    def predict_open_loop(self, dt):
        """Propagation sans mesure, quand la balle n'est pas detectee. La
        covariance croit, ce qui est le comportement voulu : l'estimation se
        degrade tant qu'aucune mesure ne la recale."""
        if not self.ready:
            return self.x.copy()
        self.x, self.P = self._step(self.x, self.P, self.u_hist[0], dt)
        x_now, P_now = self.x.copy(), self.P.copy()
        for u_k in self.u_hist[1:]:
            x_now, P_now = self._step(x_now, P_now, u_k, dt)
        self.x_now = x_now
        return x_now


def build_estimator(cfg, dt_nominal):
    """Construit l'EKF a partir d'une configuration de banc.

    Contrairement a estimator_v3.build_estimator, UNWRAP_PSI ne declenche
    aucun refus : c'est tout l'objet de ce fichier. A21 est ici le
    coefficient du rappel SINUSOIDAL, pas de sa linearisation -- les deux
    ont la meme valeur numerique (-omega_n^2), seule leur place dans le
    modele change.

    SIGMA_PSI_DEG doit correspondre au cerceau utilise : 0.4 deg mesure sur
    l'outer, 0.035 deg sur l'inner. Une variance de mesure surestimee d'un
    facteur 100 fait que le filtre ne croit plus la camera et ajoute du
    retard la ou la mesure est justement la meilleure.
    """
    n_delay = int(round(cfg.CAM_LATENCY_S / dt_nominal))
    ekf = ExtendedKF(
        a21=cfg.A21, a22=cfg.A22, b2=cfg.B2,
        sigma_psi=math.radians(cfg.SIGMA_PSI_DEG),
        q_accel=cfg.Q_ACCEL,
        n_delay=n_delay)
    return ekf, n_delay
