#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# --- Parameters (edit these) ---
L = 10.0            # channel length
H = 1.0             # channel height
N = 6               # number of baffles
x0 = 1.0            # inlet buffer before first baffle
Lb = 1.0            # baffle streamwise length support
g_min = 0.15        # minimum gap to opposite wall (validity constraint)

# Example baffle heights (must satisfy 0 < h_k < H - g_min)
h = np.array([0.45, 0.40, 0.55, 0.35, 0.50, 0.42]) * H
h = np.clip(h, 1e-3, H - g_min)

dx = (L - 2*x0)/(N-1)
xk = x0 + np.arange(N)*dx

def f_cos(x, xk, Lb):
    """Cosine bump on |x-xk|<=Lb/2, else 0."""
    s = (x - xk) / (Lb/2)
    out = np.zeros_like(x)
    m = np.abs(s) <= 1.0
    out[m] = 0.5*(1 + np.cos(np.pi*s[m]))
    return out

fig, ax = plt.subplots(figsize=(10, 2.6))
ax.set_aspect('equal', adjustable='box')

# Channel outline
ax.plot([0, L, L, 0, 0], [0, 0, H, H, 0], linewidth=2)

xx = np.linspace(0, L, 2000)
for k in range(N):
    bump = h[k] * f_cos(xx, xk[k], Lb)
    m = bump > 0
    if (k % 2) == 0:  # bottom baffle (k=0,2,4,...)
        ax.fill_between(xx[m], 0, bump[m], alpha=0.25)
        ax.plot(xx[m], bump[m], linewidth=1.5)
        ax.annotate(rf"$h_{{{k+1}}}$", xy=(xk[k], h[k]*0.5), ha="left", va="center", fontsize=10)
    else:             # top baffle
        ax.fill_between(xx[m], H - bump[m], H, alpha=0.25)
        ax.plot(xx[m], H - bump[m], linewidth=1.5)
        ax.annotate(rf"$h_{{{k+1}}}$", xy=(xk[k], H - h[k]*0.5), ha="left", va="center", fontsize=10)

    ax.annotate(rf"$x_{{{k+1}}}$", xy=(xk[k], H*1.03), ha="center", va="bottom", fontsize=10)

# Global annotations
ax.annotate(r"$\\Gamma_{\\mathrm{in}}$", xy=(0, H/2), xytext=(-0.2, H/2),
            ha="right", va="center", fontsize=11)
ax.annotate(r"$\\Gamma_{\\mathrm{out}}$", xy=(L, H/2), xytext=(L+0.2, H/2),
            ha="left", va="center", fontsize=11)

ax.text(L/2, -0.05, r"$L$", ha="center", va="top", fontsize=11)
ax.text(-0.1, H/2, r"$H$", ha="right", va="center", fontsize=11)

# Lb marker at first baffle
ax.annotate("", xy=(xk[0]-Lb/2, -0.02), xytext=(xk[0]+Lb/2, -0.02),
            arrowprops=dict(arrowstyle="<->", linewidth=1.2))
ax.text(xk[0], -0.05, r"$L_b$", ha="center", va="top", fontsize=11)

# x0 marker
ax.annotate("", xy=(0, -0.12), xytext=(x0, -0.12),
            arrowprops=dict(arrowstyle="<->", linewidth=1.2))
ax.text(x0/2, -0.15, r"$x_0$", ha="center", va="top", fontsize=11)

ax.set_xlim(-0.6, L+0.6)
ax.set_ylim(-0.25, H+0.25)
ax.axis("off")
ax.set_title("2D Two-Stream Mixer with Alternating Cosine Baffles (Parameterization Sketch)")

plt.savefig("mixer_parameterization_sketch.png", dpi=200, bbox_inches="tight")
plt.show()
