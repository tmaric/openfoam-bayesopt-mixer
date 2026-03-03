import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


def shape_function(x, xk, Lb):
    s = (x - xk) / (Lb / 2.0)
    out = 0.5 * (1.0 + np.cos(np.pi * s))
    out[np.abs(x - xk) > (Lb / 2.0)] = 0.0
    return out


def draw_mixer_sketch(output_png):
    # Geometric parameters used for illustration only.
    L = 12.0
    H = 2.0
    N = 6
    x0 = 1.5
    Lb = 1.2
    h = np.array([0.45, 0.65, 0.55, 0.75, 0.50, 0.70])

    delta_x = (L - 2.0 * x0) / (N - 1)
    xk = np.array([x0 + i * delta_x for i in range(N)])

    fig, ax = plt.subplots(figsize=(12, 4), dpi=150)

    # Channel domain boundary.
    ax.add_patch(Rectangle((0, 0), L, H, fill=False, linewidth=2.0, edgecolor='black'))

    # Alternating baffles with cosine-tip parameterization.
    for i in range(N):
        x_left = xk[i] - Lb / 2.0
        x_right = xk[i] + Lb / 2.0
        xb = np.linspace(x_left, x_right, 300)
        fk = shape_function(xb, xk[i], Lb)

        if (i + 1) % 2 == 1:
            yb = h[i] * fk
            ax.fill_between(xb, 0.0, yb, color='tab:blue', alpha=0.30, linewidth=0)
            ax.plot(xb, yb, color='tab:blue', linewidth=2.0)
        else:
            yb = H - h[i] * fk
            ax.fill_between(xb, yb, H, color='tab:orange', alpha=0.30, linewidth=0)
            ax.plot(xb, yb, color='tab:orange', linewidth=2.0)

        ax.axvline(xk[i], ymin=0.0, ymax=1.0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax.text(xk[i], H + 0.12, f"x_{i+1}", ha='center', va='bottom', fontsize=9)

    # Inlet scalar split marker.
    ax.plot([0, 0], [H / 2.0, H], color='tab:red', linewidth=3)
    ax.plot([0, 0], [0, H / 2.0], color='tab:green', linewidth=3)
    ax.text(-0.15, 0.75 * H, 'c=1', ha='right', va='center', fontsize=9, color='tab:red')
    ax.text(-0.15, 0.25 * H, 'c=0', ha='right', va='center', fontsize=9, color='tab:green')

    # Core annotations.
    ax.annotate('Inlet', xy=(0, H / 2.0), xytext=(-0.9, H / 2.0), arrowprops=dict(arrowstyle='->', lw=1.0), va='center')
    ax.annotate('Outlet', xy=(L, H / 2.0), xytext=(L + 0.8, H / 2.0), arrowprops=dict(arrowstyle='->', lw=1.0), va='center')

    ax.annotate('L', xy=(0, -0.18), xytext=(L, -0.18), arrowprops=dict(arrowstyle='<->', lw=1.0), ha='center', va='top')
    ax.annotate('H', xy=(-0.3, 0), xytext=(-0.3, H), arrowprops=dict(arrowstyle='<->', lw=1.0), ha='right', va='center')

    ax.text(
        L * 0.5,
        H + 0.35,
        r'$\theta=(h_1,\ldots,h_N),\;x_k=x_0+(k-1)\Delta x,\;\Delta x=(L-2x_0)/(N-1)$',
        ha='center',
        va='bottom',
        fontsize=10,
    )
    ax.text(
        L * 0.5,
        -0.42,
        'f_k(x)=0.5*(1+cos(pi*(x-x_k)/(L_b/2))), |x-x_k| <= L_b/2',
        ha='center',
        va='top',
        fontsize=10,
    )

    ax.set_xlim(-1.2, L + 1.2)
    ax.set_ylim(-0.6, H + 0.6)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('2D Laminar Two-Stream Passive Mixer: Geometry and Parameterization Sketch')

    fig.tight_layout()
    fig.savefig(output_png, dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    output = os.path.join(os.path.dirname(__file__), 'passive_mixer_parameterization_sketch.png')
    draw_mixer_sketch(output)
