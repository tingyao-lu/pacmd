# ------------------------------------------------------
#
# ACMD code with normal mode transformation and normal mode thermostatting
# This is to reproduce Hone, Rossky, and Voth's JCP paper in 2006
#
# ------------------------------------------------------

import numpy as np
from func import asym_pot,normal_mode
from eom import force

P = 32 # number of beads
M = 3 # NHC no.

#--------------------------------------
# index: 
# j, bead index; 1 to P
# k, normal mode index; 0 to P-1
# s, NHC index; 1 to M
#--------------------------------------

# parameters
nsteps = 10000
kb=1
hbar=1
m=1


#beta = 1/(kb*T)
beta = 8
beta_P = beta/P
n = np.arange(P-1)
omega_P = P/(beta*hbar)
omega_n = 2*omega_P*np.sin(n*np.pi/P) # internal mode frequencies
Omega_n = 
mn = m*omega_n**2/Omega_n**2 # internal mode mass
gamma = 0.4 # adiabatic separation

# beads
q = np.zeros((P,))
p = np.zeros((P,))
# centroid
pc = 
qc = np.sum(q)/P
# normal modes for beads
qj = 
C =  normal_mode(P) # normal mode transformation matrix Cjk
Q =  C.T @ q # normal mode pos
P =  # internal mode mom


rNHC = 
pNHC = 
mNHC = 


Zc = np.sum(np.exp(-beta_P * Ham)) if qc == np.sum(q)/P else None

pmf = np.ln(Zc)/beta


# initialization
# MB sampling of initial velocity
# 

for i in range(nsteps):
    qj = qc + np.sum(Cjn*Qn)
    Fj = force(asym_pot(qj, c, g), qj) #bead force

    ham(p,m,m,P,lambda,a,pot,U,hbar,beta)

    Fc = np.sum(Fj)/P # centroid force
    Fn = -m*omega_n**2*Qn + np.sum(Cjn*Fj)/P # normal mode force 

    # propagation
    




