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

def pot(q):
    return 0.25*q**4

def forcequad(q):
    return -q**3

def forcespr(omega, m,P):
    return -omega**2*m*P


def vvint(q,v,f, m, dt, force):
    q_new = q + dt*v+dt**2*f/(2*m)
    f_new = force(q_new)
    v_new = v + dt*(f+f_new)/(2*m)
    p_new = m*v_new
    return q_new, f_new, v_new, p_new
    

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
    