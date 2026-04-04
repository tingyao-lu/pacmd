import numpy as np
import matplotlib.pyplot as plt

data = np.load("pacmd_output.npz")
time        = data["time"]
q_centroid  = data["q_centroid"]
energy      = data["energy"]
corr_xx     = data["corr_xx"]
dt          = float(data["dt"])
n_equil     = int(data["n_equil"])

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Panel 1 — centroid position
axes[0].plot(time, q_centroid)
axes[0].set_xlabel("Time")
axes[0].set_ylabel(r"$q_c(t)$")
axes[0].set_title("Centroid position")

# Panel 2 — energy conservation
axes[1].plot(time, energy)
axes[1].set_xlabel("Time")
axes[1].set_ylabel("Energy")
axes[1].set_title("Energy conservation")

# Panel 3 — position autocorrelation (post-equilibration)
t_corr = np.arange(len(corr_xx)) * dt
axes[2].plot(t_corr, corr_xx)
axes[2].set_xlabel("Lag time")
axes[2].set_ylabel(r"$C_{xx}$")
axes[2].set_title("Centroid position autocorrelation")

fig.tight_layout()
plt.savefig("plot_pacmd_results.png", dpi=150)
plt.show()
