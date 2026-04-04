import numpy as np
import matplotlib.pyplot as plt
from func import normal_mode, forcequad, pot, Cxx

# ------------------------------------------------------
# PACMD core (no thermostat yet)
# Assumptions:
#   1. normal_mode(P) returns an orthonormal matrix C of shape (P, P)
#   2. bead <-> normal-mode transforms are:
#          Q = C.T @ q_beads
#          q_beads = C @ Q
#   3. forcequad(q) returns the PHYSICAL FORCE -dV/dq, vectorized over q
#   4. pot(q) returns the physical potential V(q), vectorized over q
# ------------------------------------------------------

# parameters
P = 32                  # number of beads = number of normal modes
nsteps = 10000
n_equil = 2000          # throw away this many steps before analysis
kB = 1.0
hbar = 1.0
m = 1.0
beta = 8.0
beta_P = beta / P
dt = 0.01
rng = np.random.default_rng(1234)

# free ring-polymer normal-mode frequencies
k = np.arange(P)  # 0 ... P-1
omega_k = (2.0 * P / (beta * hbar)) * np.sin(np.pi * k / P)
omega_k[0] = 0.0

# PACMD target internal-mode frequency
Omega = P ** (P / (P - 1)) / (beta * hbar)

# PACMD fictitious masses
m_tilde = np.empty(P)
m_tilde[0] = m
m_tilde[1:] = m * (omega_k[1:] ** 2) / (Omega ** 2)

# normal-mode matrix
C = normal_mode(P)
assert C.shape == (P, P)
assert np.allclose(C.T @ C, np.eye(P), atol=1e-12), "C must be orthonormal"

# storage
Q_nm = np.zeros((nsteps, P))      # normal-mode positions
P_nm = np.zeros((nsteps, P))      # normal-mode momenta
q_beads = np.zeros((nsteps, P))   # bead positions
energy = np.zeros(nsteps)

# ---------- transforms ----------
def beads_to_nm(qb):
    return C.T @ qb

def nm_to_beads(Q):
    return C @ Q

# ---------- forces ----------
def physical_force_nm(Q):
    qb = nm_to_beads(Q)
    F_beads = forcequad(qb)
    return beads_to_nm(F_beads)

def spring_force_nm(Q):
    return -m * (omega_k ** 2) * Q

def total_force_nm(Q):
    return physical_force_nm(Q) + spring_force_nm(Q)

# ---------- energy ----------
def pacmd_energy(Q, Pnm):
    qb = nm_to_beads(Q)
    kinetic = np.sum(Pnm**2 / (2.0 * m_tilde))
    spring = 0.5 * np.sum(m * (omega_k**2) * Q**2)
    potential = np.sum(pot(qb))
    return kinetic + spring + potential

# ------------------------------------------------------
# initialization
# ------------------------------------------------------

# simplest initial coordinates: all beads at 0
q_beads[0] = np.zeros(P)
Q_nm[0] = beads_to_nm(q_beads[0])
P_nm[0] = rng.normal(loc=0.0, scale=np.sqrt(m_tilde / beta_P), size=P)
P_nm[0, 0] = 0.0 # centroid

# initial force and energy
F_nm = total_force_nm(Q_nm[0])
energy[0] = pacmd_energy(Q_nm[0], P_nm[0])

# ------------------------------------------------------
# propagate the beads in normal-mode space using velocity Verlet
# ------------------------------------------------------
for i in range(nsteps - 1):
    # position step
    Q_next = Q_nm[i] + (P_nm[i] / m_tilde) * dt + 0.5 * (F_nm / m_tilde) * dt**2    

    # new force
    F_next = total_force_nm(Q_next)

    # momentum step
    P_next = P_nm[i] + 0.5 * (F_nm + F_next) * dt

    # store
    Q_nm[i + 1] = Q_next
    P_nm[i + 1] = P_next
    q_beads[i + 1] = nm_to_beads(Q_next)
    energy[i + 1] = pacmd_energy(Q_next, P_next)

    F_nm = F_next

# ------------------------------------------------------
# observables
# ------------------------------------------------------
q_centroid = q_beads.mean(axis=1)
corr_xx = Cxx(q_centroid[n_equil:], dt)

time = np.arange(nsteps) * dt

plt.figure()
plt.plot(time, q_centroid, label="centroid position")
plt.xlabel("time")
plt.ylabel(r"$q_c(t)$")
plt.title("Centroid position")
plt.legend()
plt.show()

plt.figure()
plt.plot(time, energy, label="PACMD energy")
plt.xlabel("time")
plt.ylabel("energy")
plt.title("Energy conservation check")
plt.legend()
plt.show()

plt.figure()
plt.plot(np.arange(len(corr_xx)) * dt, corr_xx, label="Cxx")
plt.xlabel("time")
plt.ylabel("correlation")
plt.title("Centroid position autocorrelation")
plt.legend()
plt.show()
