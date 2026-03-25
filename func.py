import numpy as np

def pot(x):
    return 0.25*x**4

def maxwell_boltzmann_speed(v, m=1.0, T=1.0, kB=1.0):
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
    coef = np.sqrt((m / (2 * np.pi * kB * T))**3)
    return 4 * np.pi * coef * v**2 * np.exp(-m * v**2 / (2 * kB * T))


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
    for m in range(1, max_m + 1):
        for j in range(P):
            angle = 2.0 * np.pi * m * j / P
            C[j, col] = np.sqrt(2.0 / P) * np.cos(angle)
            C[j, col + 1] = np.sqrt(2.0 / P) * np.sin(angle)
        col += 2

    # extra alternating mode for even P
    if P % 2 == 0:
        for j in range(P):
            C[j, P - 1] = (-1)**j / np.sqrt(P)

    return C


def ham(p,m,m,P,lambda,a,pot,U,): # normal mode Hamiltonian
    '''
    vector: p_k, m_k, lambda, a, pot
    matrix: U
    parameter: m, hbar, beta
    returns: a matrix
    '''
    p1 = np.square(p)/(2*m)
    p2 = (m*P/(2*hbar**2*beta**2))*np.multiply(lambda, np.square(a))
    p3 = pot(np.sqrt(P)* U@a)/P
    return p1+p2+p3

def mbdis():
    return


