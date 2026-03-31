# ------------------------------------------------------
#
# PACMD code with normal mode transformation and normal mode thermostatting
# This is to reproduce Hone, Rossky, and Voth's JCP paper in 2006
#
# ------------------------------------------------------

import numpy as np
from func import forcespr, forcequad, normal_mode, pot, vvint
from openmm import NoseHooverChain
import matplotlib.pyplot as plt

# parameters
P = 32 # beads no.
M = 3 # NHC no.
nsteps = 10000
kb=1
hbar=1
m=1
beta = 8  #beta = 1/(kb*T)
T=1/(kb*beta) # 0.125 a.u.
beta_P = beta/P
gamma = 0.4 # adiabatic separation
dt=0.01

k = np.arange(0,P+1) # normal mode index, 0 to P
omega_k = 2*P/(beta*hbar)*np.sin(k*np.pi/P)  # internal mode frequencies [P, ], eq23
omega_k[0] = 0 # set k=0 mode frequency to zero to avoid singularity
Omega = P**(P/(P-1))/(beta*hbar) # scalar

m_k = m * omega_k**2 / Omega**2 # [P+1, ]
m_k[0] = m # set k=0 mode mass to physical mass

# centroid and beads positions and momenta
q, p = np.zeros((nsteps,P+1)), np.zeros((nsteps,P+1)) # [nsteps,P+1]

# normal modes for beads
C =  normal_mode(P) # normal mode transformation matrix Cjk [P,P]
# Vectorized transform for all timesteps: (nsteps, P) @ (P, P) -> (nsteps, P)
Q = q[:,1:] @ C   # [nsteps, P]
Pnm = p[:,1:] @ C # [nsteps, P]

# forces
f = np.zeros((nsteps,P+1)) #[nsteps,P+1]

#----------------
# initialization
#----------------
q_init = np.zeros(P) #[P,]
p_init = m_k[1:]*np.random.normal(0,1/np.sqrt(beta_P*m_k[1:]),P) #[P,]

q_init[0] = 0 # set centroid position to zero
p_init[0] = 0 # set centroid momentum to zero

#beads
p[0,1:] = p_init
q[0,1:] = q_init
Pnm[0,:] = C.T @ p[0,1:]
Q[0,:] = C.T @ q[0,1:]

#centroid
p[0,0]=np.sum(p[0,1:])/P
q[0,0]=np.sum(q[0,1:])/P

#force = f(ext)+f(spr)
f[0,0] = forcequad(q[0,0]) + forcespr(omega_k[0], m_k[0], p[0,0])
f[0,1:] = forcequad(Q[0,:]) + forcespr(omega_k[1:], m_k[1:], Pnm[0,:])


#----------------
# propagation
#----------------
i = 0
for i in range(nsteps-1):
    #beads
    Q[i+1,:], f[i+1,1:], _, Pnm[i+1,:] = vvint(
        Q[i,:], Pnm[i,:]/m_k[1:],f[i,1:],m_k[1:],dt,forcequad(Q[i,:]) + forcespr(omega_k[1:], m_k[1:], Pnm[i,:]))

    # Keep centroid mode separate from internal normal modes.

    p[i+1,1:] = C @ Pnm[i+1,:]
    q[i+1,1:] = C @ Q[i+1,:]

    if i == 0:
        print('i = 1')
        #print(Q[i+1,:])

    #centroid
    q[i+1,0] = np.sum(q[i+1,1:])/P
    p[i+1,0] = np.sum(p[i+1,1:])/P
    f[i+1,0] = forcequad(q[i+1,0])


    Pnm[i+1,:] = C.T @ p[i+1,1:]
    Q[i+1,:] = C.T @ q[i+1,1:]

    i += 1

#print(qc)
steplist = np.arange(nsteps)
#plt.plot(steplist, pot(steplist))
plt.plot(steplist, q[:,0])
plt.xlabel('Time step')
plt.ylabel('Centroid position')
plt.title('Centroid Position vs Time')
plt.show()







