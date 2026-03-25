
def pot(x):
    return 0.25*x**4

def ImgTimeCorrFunc():
    return 

def RealTimeCorrFunc():
    return 

def KuboCorrFunc():
    return 


def ham(p,m,m,P,lambda,a,pot,U,hbar,beta): # normal mode Hamiltonian
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


