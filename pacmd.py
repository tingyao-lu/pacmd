import numpy as np
import matplotlib.pyplot as plt
from func import (normal_mode, forcequad, pot, Cxx, beads_to_nm, nm_to_beads,
                  physical_force_nm, spring_force_nm, total_force_nm,
                  pacmd_energy, vvint, NHC)


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

# NHC thermostat for internal modes (k=1..P-1)
# Target temperature for each mode: kT = 1/beta_P
# M_chain = 3                                      # chain length
# tau_nhc = 1.0 / Omega                            # coupling timescale ~ 1/Omega
# kT_target = np.full(P, kB / beta_P)              # same target kT for all modes
# nhc = NHC(M_chain, P, m_tilde, kT_target, tau_nhc, n_mts=1, n_ys=3)


# ------------------------------------------------------
# initialization
# ------------------------------------------------------

# simplest initial coordinates: all beads at 0
q_beads[0] = np.zeros(P)
Q_nm[0] = beads_to_nm(P, q_beads[0])
P_nm[0] = rng.normal(loc=0.0, scale=np.sqrt(m_tilde / beta_P), size=P)
P_nm[0, 0] = 0.0 # centroid

# initial force and energy
F_nm = total_force_nm(P, m, omega_k, Q_nm[0])
energy[0] = pacmd_energy(P, m, m_tilde, omega_k, Q_nm[0], P_nm[0])

# ------------------------------------------------------
# propagate the beads in normal-mode space using velocity Verlet
# ------------------------------------------------------
for i in range(nsteps - 1):
    # --- NHC half-step: thermostat internal momenta ---
    P_cur = P_nm[i].copy()
    #P_cur[1:] = nhc.step(P_cur[1:], dt / 2)

    # --- velocity Verlet position step ---
    Q_next = Q_nm[i] + (P_cur / m_tilde) * dt + 0.5 * (F_nm / m_tilde) * dt**2

    # new force
    F_next = total_force_nm(P, m, omega_k, Q_next)

    # --- velocity Verlet momentum step ---
    P_next = P_cur + 0.5 * (F_nm + F_next) * dt

    # --- NHC half-step: thermostat internal momenta ---
    #P_next[1:] = nhc.step(P_next[1:], dt / 2)

    # store
    Q_nm[i + 1] = Q_next
    P_nm[i + 1] = P_next
    q_beads[i + 1] = nm_to_beads(P, Q_next)
    energy[i + 1] = pacmd_energy(P, m, m_tilde, omega_k, Q_next, P_next)

    F_nm = F_next

# ------------------------------------------------------
# observables
# ------------------------------------------------------
q_centroid = q_beads.mean(axis=1)
corr_xx = Cxx(q_centroid[n_equil:], dt)

time = np.arange(nsteps) * dt

np.savez("pacmd_output.npz",
         time=time, q_centroid=q_centroid, energy=energy,
         corr_xx=corr_xx, dt=dt, n_equil=n_equil)
