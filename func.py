import numpy as np
hbar = 1

def mbdis(v, m=1.0, T=1.0, kB=1.0):
    """Maxwell-Boltzmann speed distribution f(v) for a gas of particles.

    Parameters
    ----------
    v : array_like or float
        Speed(s) at which to evaluate the distribution.
    m : float
        Particle mass.
    T : float
        Temperature.
    kB : float
        Boltzmann constant.

    Returns
    -------
    f : ndarray or float
        Probability density f(v).
    """
    v = np.asarray(v)
    coef = np.sqrt((m / (2 * np.pi * kB * T))**1/2)
    return 4 * np.pi * coef * v**2 * np.exp(-m * v**2 / (2 * kB * T))

def Zk(k, P, betaP):
    np.exp()
    return 1/(2*np.pi*hbar)

# def Ham(p, q, m, P, lambda_, a, pot, U, hbar, beta):
#     term1 = np.square(p)/(2*m)
#     term2 = 
#     term3 = 
#     return term1 + term2 + term3

def ImgTimeCorrFunc():
    return 

def RealTimeCorrFunc():
    return 

def KuboCorrFunc():
    return 

#--------------------------
# normal mode transforms
#--------------------------

def normal_mode(P):
    """
    Build the real-valued ring-polymer normal-mode transform matrix C.

    Q_k = sum_j C[j,k] q_j
    q_j = sum_k C[j,k] Q_k

    Returns
    -------
    C : ndarray, shape (P, P)
        Orthogonal normal-mode transform matrix
    """
    C = np.zeros((P, P))

    # centroid mode
    for j in range(P):
        C[j, 0] = 1.0 / np.sqrt(P)

    col = 1
    max_m = (P - 1) // 2

    # cosine/sine pairs
    for n in range(1, max_m + 1):
        for j in range(P):
            angle = 2.0 * np.pi * n * j / P
            C[j, col] = np.sqrt(2.0 / P) * np.cos(angle)
            C[j, col + 1] = np.sqrt(2.0 / P) * np.sin(angle)
        col += 2

    # extra alternating mode for even P
    if P % 2 == 0:
        for j in range(P):
            C[j, P - 1] = (-1)**j / np.sqrt(P)

    return C

def beads_to_nm(P,qb):
    C = normal_mode(P)
    return C.T @ qb

def nm_to_beads(P,Q):
    C = normal_mode(P)
    return C @ Q

#--------------------------
# forces
#--------------------------

def physical_force_nm(P, Q):
    qb = nm_to_beads(P, Q)
    F_beads = forcequad(qb)
    return beads_to_nm(P, F_beads)

def spring_force_nm(m, omega_k, Q):
    return -m * (omega_k ** 2) * Q

def total_force_nm(P, m, omega_k, Q):
    return physical_force_nm(P, Q) + spring_force_nm(m, omega_k, Q)


#--------------------------
# energies
#--------------------------

def pacmd_energy(P, m, m_tilde, omega_k, Q, Pnm):
    qb = nm_to_beads(P, Q)
    kinetic = np.sum(Pnm**2 / (2.0 * m_tilde))
    spring = 0.5 * np.sum(m * (omega_k**2) * Q**2)
    potential = np.sum(pot(qb))
    return kinetic + spring + potential


def pot(q):
    return 0.25*q**4

def forcequad(q):
    return -q**3

def forcespr(omega, m,P):
    return -omega**2*m*P


def Cxx(q, dt):
    q = q-np.mean(q, axis=0) # remove mean
    nsteps = len(q)
    C = np.zeros(nsteps)
    for k in range(nsteps):
        C[k] = np.mean(q[:nsteps-k] * q[k:])
    return C/C[0]

def vvint(q,v,f, m, dt, force):
    q_new = q + dt*v+dt**2*f/(2*m)
    f_new = force(q_new)
    v_new = v + dt*(f+f_new)/(2*m)
    p_new = m*v_new
    return q_new, f_new, v_new, p_new

#-------------------
# Nose-Hoover chain thermostat
# Ref: Martyna, Tuckerman, Tobias & Klein, Mol. Phys. 87, 1117 (1996)
#-------------------

# Yoshida-Suzuki weights for the symmetric decomposition
_YS_WEIGHTS = {
    1: np.array([1.0]),
    3: np.array([1.351207191959657, -1.702414383919315, 1.351207191959657]),
    5: np.array([0.41449077179437573, 0.41449077179437573,
                 -0.65796308717750293, 0.41449077179437573,
                 0.41449077179437573]),
}


class NHC:
    """Nose-Hoover Chain thermostat for PACMD internal normal modes.

    Each internal mode k gets its own chain of M thermostat variables.
    The centroid (k=0) is NOT thermostatted.

    Parameters
    ----------
    M : int
        Chain length (number of thermostat DOFs per mode).
    P : int
        Number of normal modes (= number of beads).
    m_tilde : ndarray, shape (P,)
        Fictitious masses for each normal mode.
    kT_target : ndarray, shape (P,)
        Target kinetic energy per mode = kB * T_target for each mode.
        For PACMD internal modes: kT_target[k] = 1 / beta_P.
    tau : float
        Thermostat coupling time scale.
    n_mts : int
        Number of multi-timestep iterations (n_c in the paper).
    n_ys : int
        Yoshida-Suzuki order (1, 3, or 5).
    """

    def __init__(self, M, P, m_tilde, kT_target, tau, n_mts=1, n_ys=3):
        self.M = M
        self.P = P
        self.n_modes = P - 1            # internal modes 1..P-1
        self.m_tilde = m_tilde[1:]      # shape (P-1,)
        self.kT = kT_target[1:]         # shape (P-1,)
        self.n_mts = n_mts
        self.n_ys = n_ys
        self.w = _YS_WEIGHTS[n_ys]

        # thermostat "masses"  Q[s] for each mode and chain link
        # Q[s=0] = kT * tau^2,  Q[s>0] = kT * tau^2  (common choice)
        # shape: (P-1, M)
        self.Q = np.outer(self.kT, np.full(M, tau**2))

        # thermostat momenta and positions — one chain per internal mode
        self.p_eta = np.zeros((self.n_modes, M))  # (P-1, M)
        self.eta   = np.zeros((self.n_modes, M))  # (P-1, M)

    def step(self, P_nm_internal, dt):
        """Apply one full NHC thermostat step (half-step before + after VV).

        Call this as a half-step: once before the Verlet position update
        and once after the Verlet momentum update — or equivalently wrap
        the Verlet step with nhc.step(p, dt/2) on each side.

        Parameters
        ----------
        P_nm_internal : ndarray, shape (P-1,)
            Momenta of internal normal modes (modified in-place).
        dt : float
            The timestep (pass dt/2 for a symmetric split).

        Returns
        -------
        P_nm_internal : ndarray, shape (P-1,)
            Thermostatted momenta (same array, modified in-place).
        """
        M = self.M
        kT = self.kT                    # (P-1,)
        Q  = self.Q                     # (P-1, M)
        p_eta = self.p_eta              # (P-1, M)
        m_t = self.m_tilde              # (P-1,)

        for _ in range(self.n_mts):
            for w_j in self.w:
                delta = w_j * dt / self.n_mts

                # current kinetic energy per mode
                KE = P_nm_internal**2 / m_t       # (P-1,)

                # "force" on the first thermostat variable
                G = KE - kT                        # (P-1,)

                # propagate chain from the deepest link inward
                # --- last link M-1 ---
                p_eta[:, M-1] += 0.25 * delta * self._G_chain(M-1, G, p_eta, Q, kT)

                # --- links s = M-2 .. 1 ---
                for s in range(M-2, 0, -1):
                    xi = p_eta[:, s+1] / Q[:, s+1]
                    p_eta[:, s] *= np.exp(-0.125 * delta * xi)
                    p_eta[:, s] += 0.25 * delta * self._G_chain(s, G, p_eta, Q, kT)
                    p_eta[:, s] *= np.exp(-0.125 * delta * xi)

                # --- link s=0 (couples to the physical momentum) ---
                if M > 1:
                    xi = p_eta[:, 1] / Q[:, 1]
                    p_eta[:, 0] *= np.exp(-0.125 * delta * xi)
                p_eta[:, 0] += 0.25 * delta * G
                if M > 1:
                    p_eta[:, 0] *= np.exp(-0.125 * delta * xi)

                # scale physical momenta
                xi0 = p_eta[:, 0] / Q[:, 0]
                scale = np.exp(-0.5 * delta * xi0)
                P_nm_internal *= scale

                # update thermostat positions
                self.eta += 0.5 * delta * p_eta / Q

                # recompute G after scaling
                KE = P_nm_internal**2 / m_t
                G = KE - kT

                # propagate chain forward (mirror of backward sweep)
                # --- link s=0 ---
                if M > 1:
                    xi = p_eta[:, 1] / Q[:, 1]
                    p_eta[:, 0] *= np.exp(-0.125 * delta * xi)
                p_eta[:, 0] += 0.25 * delta * G
                if M > 1:
                    p_eta[:, 0] *= np.exp(-0.125 * delta * xi)

                # --- links s=1 .. M-2 ---
                for s in range(1, M-1):
                    xi = p_eta[:, s+1] / Q[:, s+1]
                    p_eta[:, s] *= np.exp(-0.125 * delta * xi)
                    p_eta[:, s] += 0.25 * delta * self._G_chain(s, G, p_eta, Q, kT)
                    p_eta[:, s] *= np.exp(-0.125 * delta * xi)

                # --- last link ---
                p_eta[:, M-1] += 0.25 * delta * self._G_chain(M-1, G, p_eta, Q, kT)

        return P_nm_internal

    @staticmethod
    def _G_chain(s, G_phys, p_eta, Q, kT):
        """Compute the 'force' on chain link s.

        s=0:  G = KE_phys - kT  (passed as G_phys).
        s>0:  G = p_eta[s-1]^2 / Q[s-1] - kT.
        """
        if s == 0:
            return G_phys
        return p_eta[:, s-1]**2 / Q[:, s-1] - kT



    

# def ham(p,m,m,P,lambda,a,pot,U,): # normal mode Hamiltonian
#     '''
#     vector: p_k, m_k, lambda, a, pot
#     matrix: U
#     parameter: m, hbar, beta
#     returns: a matrix
#     '''
#     p1 = np.square(p)/(2*m)
#     p2 = (m*P/(2*hbar**2*beta**2))*np.multiply(lambda, np.square(a))
#     p3 = pot(np.sqrt(P)* U@a)/P
#     return p1+p2+p3

# def hamHarmSprNM(P, lambda, Q, m, hbar, beta):
#     return m*P*lambda*Q/(2*hbar**2*beta**2)

# def F_HarmSprNM(P, lambda, Qold, Qnew, dq, m, hbar, beta):
#     return -(hamHarmSprNM(P, lambda, Qnew, m, hbar, beta)-hamHarmSprNM(P, lambda, Qold, m, hbar, beta))/2*dq    

# def pot(q):
#     return 0.25*q**4

# def F_physNM(P, Qold, Qnew, dq, C):
#     return -(pot(np.sqrt(P)*C@Qold)-pot(np.sqrt(P)*C@Qnew))/2*dq

# def nhc_eom():
#     return

# def force(potential, q):
#     return -np.gradient(potential(q))

# def forcespr(omega, P):
#     return -omega**2*P


# def vvint(q,v,f, m, dt):
#     q_new = q + dt*v+dt**2*f/(2*m)
#     f_new = f(q_new)
#     v_new = v + dt*(f+f_new)/(2*m)
#     p_new = m*v_new
#     return q_new, f_new, v_new, p_new
    