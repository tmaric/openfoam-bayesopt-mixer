#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# --- SAR lamination ladder sketch parameters (edit as needed) ---
H = 1.0               # channel height (normalized)
L0 = 2.0              # inlet/outlet buffer length
N = 5                 # number of unit cells
L_cell = 4.0          # unit cell length
L = 2*L0 + N*L_cell   # total length

# Geometry parameters (theta)
w_s = 0.38*H          # nominal subchannel width after split
t_s = 0.10*H          # splitter thickness
L_s = 1.4             # split length within cell
L_m = 1.0             # merge length within cell
delta = 0.08*H        # top deflector vertical bias (conceptual)
r = 0.03*H            # fillet radius (annotated, not explicitly drawn)

t_m = 0.05*H          # merge splitter thickness (kept fixed for sketch)
L_c = L_cell - L_s - L_m
assert L_c > 0, "Choose L_cell > L_s + L_m"

# Derived deflector height (amount intruding from walls)
h_d = 0.5*H - w_s
assert h_d > 0, "Need w_s < H/2"

def g_env(xi, Lc):
    """Cosine envelope on [0,Lc]."""
    return 0.5*(1 - np.cos(2*np.pi*xi/Lc))

# Build plot
fig, ax = plt.subplots(figsize=(11, 2.8))
ax.set_aspect('equal', adjustable='box')

# Channel outline
ax.plot([0, L, L, 0, 0], [0, 0, H, H, 0], linewidth=2)

# Draw unit cells with split/shuffle/merge obstacles
for k in range(N):
    xk = L0 + k*L_cell

    # Splitter rectangle (split section)
    xs0, xs1 = xk, xk + L_s
    y0s, y1s = (H - t_s)/2, (H + t_s)/2
    ax.fill([xs0, xs1, xs1, xs0], [y0s, y0s, y1s, y1s], alpha=0.25)
    ax.plot([xs0, xs1, xs1, xs0, xs0], [y0s, y0s, y1s, y1s, y0s], linewidth=1.2)

    # Shuffle deflectors (cosine bumps from bottom and top)
    xc0, xc1 = xk + L_s, xk + L_s + L_c
    xx = np.linspace(xc0, xc1, 400)
    xi = xx - xc0
    bump = h_d * g_env(xi, L_c)

    # bottom deflector
    ax.fill_between(xx, 0, bump, alpha=0.25)
    ax.plot(xx, bump, linewidth=1.2)

    # top deflector with conceptual bias delta (clipped)
    top_curve = H - np.clip(bump + delta, 0, H)
    ax.fill_between(xx, top_curve, H, alpha=0.25)
    ax.plot(xx, top_curve, linewidth=1.2)

    # Merge splitter rectangle (merge section)
    xm0, xm1 = xk + L_s + L_c, xk + L_cell
    y0m, y1m = (H - t_m)/2, (H + t_m)/2
    ax.fill([xm0, xm1, xm1, xm0], [y0m, y0m, y1m, y1m], alpha=0.25)
    ax.plot([xm0, xm1, xm1, xm0, xm0], [y0m, y0m, y1m, y1m, y0m], linewidth=1.2)

    # Cell delimiter (light)
    ax.plot([xk, xk], [0, H], linewidth=0.6, alpha=0.3)
    ax.plot([xk+L_cell, xk+L_cell], [0, H], linewidth=0.6, alpha=0.3)

# Annotations
ax.annotate(r"$\\Gamma_{\\mathrm{in}}$", xy=(0, H/2), xytext=(-0.35, H/2),
            ha="right", va="center", fontsize=11)
ax.annotate(r"$\\Gamma_{\\mathrm{out}}$", xy=(L, H/2), xytext=(L+0.35, H/2),
            ha="left", va="center", fontsize=11)

ax.text(L/2, -0.08, r"$L$", ha="center", va="top", fontsize=11)
ax.text(-0.1, H/2, r"$H$", ha="right", va="center", fontsize=11)

# Parameter callouts
ax.text(L0 + 0.3, (H+t_s)/2 + 0.07, r"$t_s$", ha="left", va="bottom", fontsize=10)
ax.text(L0 + 0.7, H*0.6, r"$w_s$", ha="left", va="center", fontsize=10)
ax.text(L0 + L_s + 0.1, h_d*0.8, r"$h_d$", ha="left", va="center", fontsize=10)
ax.text(L0 + L_s + 0.1, H - h_d*0.8, r"$\\delta$", ha="left", va="center", fontsize=10)
ax.text(L0 + 0.2, H*1.03, r"$L_0$", ha="left", va="bottom", fontsize=11)
ax.text(L0 + 0.2, H*0.92, r"$N,\\,L_{\\mathrm{cell}}$", ha="left", va="top", fontsize=10)

ax.set_xlim(-0.8, L+0.8)
ax.set_ylim(-0.25, H+0.25)
ax.axis("off")
ax.set_title("SAR Lamination Ladder Mixer (2D Unit-Cell Sketch)")

png_path = "sar_lamination_ladder_mixer_sketch.png"
plt.savefig(png_path, dpi=200, bbox_inches="tight")
plt.show()
