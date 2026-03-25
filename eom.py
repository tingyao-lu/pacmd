
def eom_x(p,m):
    return p/m

def eom_p(E,m):
    return -E/m

def F_phys(V,x):
    return -np.gradient(V,x)

def F_harNM(m, omega, q):
    return -m*omega**2*q

def F_nhc():
    return 

