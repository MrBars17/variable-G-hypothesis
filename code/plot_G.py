"""
Plot the phenomenological G(t) curve:

  •   G = 0        at t = 0
  •   G = 0.85 G0  by t = 200 s  (BBN plateau)
  •   Logistic rise to 0.98 G0   by t = 380 kyr (recombination)
  •   Logistic rise to 1.00 G0   by t = 13.8 Gyr (today)
  •   Slow drift to 1.02 G0      by t = 27.6 Gyr (Big‑Hole pre‑phase)

Outputs:  figures/G_evolution.png  (log‑time axis, today marked)
"""

import numpy as np
import matplotlib.pyplot as plt

# --- timeline in seconds ---
t_BBN_end = 200.0                      # 200 s
t_rec   = 3.8e5 * 3.154e7              # 380 kyr
t_today = 13.8e9 * 3.154e7             # 13.8 Gyr
t_future = (13.8e9 + 13.8e9) * 3.154e7 # +13.8 Gyr  (27.6 Gyr total)

# --- plateau levels (fraction of G0) ---
k1, k2, k3, k4 = 0.85, 0.98, 1.00, 1.02

# --- timescales for smooth exponentials ---
tauA = 80.0                            # 80 s
tauB = (t_rec - t_BBN_end) / 3.0
tauC = (t_today - t_rec) / 3.0
tauD = 5.0e9 * 3.154e7                 # 5 Gyr

def G_frac(t):
    """Return G(t)/G0."""
    if t < t_BBN_end:
        return k1 * (1 - np.exp(-t / tauA))
    elif t < t_rec:
        return k1 + (k2 - k1) * (1 - np.exp(-(t - t_BBN_end) / tauB))
    elif t < t_today:
        return k2 + (k3 - k2) * (1 - np.exp(-(t - t_rec) / tauC))
    else:
        return k3 + (k4 - k3) * (1 - np.exp(-(t - t_today) / tauD))

G_frac_vec = np.vectorize(G_frac)

# --- time grid ---
t = np.logspace(-4, np.log10(t_future), 6000)  # seconds
Gvals = G_frac_vec(t)

# --- plotting ---
plt.figure(figsize=(7, 4))
plt.plot(t / 3.154e7 / 1e9, Gvals, lw=2, color="tab:orange")
plt.scatter([13.8], [1.0], color="red", s=30, zorder=5, label="today")
plt.xscale("log")
plt.ylim(0, 1.04)
plt.xlabel("Cosmic time [Gyr]")
plt.ylabel("$G(t) / G_0$")
plt.title("Hypothetical evolution of Newton's constant")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig("figures/G_evolution.png", dpi=300)
print("Saved figures/G_evolution.png")