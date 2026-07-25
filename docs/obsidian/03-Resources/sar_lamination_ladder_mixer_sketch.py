#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

# --- SAR lamination ladder sketch parameters (edit as needed) ---
H = 1.0               # channel height (normalized)
L0 = 2.0              # inlet/outlet buffer length
N = 5                 # number of unit cells
L_cell = 4.0          # unit cell length
L = 2 * L0 + N * L_cell

# Geometry parameters (theta)
w_s = 0.38 * H        # nominal subchannel width after split
t_s = 0.10 * H        # splitter thickness
L_s = 1.4             # split length within cell
L_m = 1.0             # merge length within cell
delta = 0.08 * H      # top deflector vertical bias (conceptual)
r = 0.03 * H          # fillet radius (annotated, not explicitly drawn)

t_m = 0.05 * H        # merge splitter thickness (kept fixed for sketch)
L_c = L_cell - L_s - L_m
assert L_c > 0, "Choose L_cell > L_s + L_m"

h_d = 0.5 * H - w_s
assert h_d > 0, "Need w_s < H/2"

Y_SCALE = 2.0
H_plot = H * Y_SCALE

FIG_BG = "#f4efe7"
PANEL_BG = "#fffdf8"
TEXT = "#15324a"
GUIDE = "#8ea6ba"
GUIDE_SOFT = "#d4dde6"
CHANNEL_EDGE = "#2a6f9e"
SPLIT_FACE = "#f6b26b"
SPLIT_EDGE = "#b45f06"
BOTTOM_FACE = "#7db2d9"
BOTTOM_EDGE = "#2f6f9f"
TOP_FACE = "#f7a1a1"
TOP_EDGE = "#c94b4b"
MERGE_FACE = "#93c47d"
MERGE_EDGE = "#3c8d2f"


def ys(values):
    return np.asarray(values) * Y_SCALE


def g_env(xi, Lc):
    return 0.5 * (1 - np.cos(2 * np.pi * xi / Lc))


def dim_h(ax, x0, x1, y, label, text_dy=0.10):
    ax.annotate(
        "",
        xy=(x1, y),
        xytext=(x0, y),
        arrowprops=dict(arrowstyle="<->", lw=1.4, color=GUIDE, shrinkA=0, shrinkB=0),
    )
    ax.plot([x0, x0], [y - 0.06, y + 0.06], color=GUIDE, lw=1.2)
    ax.plot([x1, x1], [y - 0.06, y + 0.06], color=GUIDE, lw=1.2)
    ax.text(
        0.5 * (x0 + x1),
        y + text_dy,
        label,
        ha="center",
        va="bottom",
        fontsize=12,
        color=TEXT,
        bbox=dict(boxstyle="round,pad=0.18", fc=PANEL_BG, ec="none"),
    )


def dim_v(ax, x, y0, y1, label, text_dx=0.13):
    ax.annotate(
        "",
        xy=(x, y1),
        xytext=(x, y0),
        arrowprops=dict(arrowstyle="<->", lw=1.4, color=GUIDE, shrinkA=0, shrinkB=0),
    )
    ax.plot([x - 0.05, x + 0.05], [y0, y0], color=GUIDE, lw=1.2)
    ax.plot([x - 0.05, x + 0.05], [y1, y1], color=GUIDE, lw=1.2)
    ax.text(
        x - text_dx,
        0.5 * (y0 + y1),
        label,
        ha="right",
        va="center",
        fontsize=12,
        color=TEXT,
        bbox=dict(boxstyle="round,pad=0.18", fc=PANEL_BG, ec="none"),
    )


def callout(ax, text, xy, xytext, ha="center", va="center", rad=0.0):
    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        ha=ha,
        va=va,
        fontsize=12,
        color=TEXT,
        bbox=dict(boxstyle="round,pad=0.28", fc=PANEL_BG, ec=GUIDE_SOFT, lw=1.0),
        arrowprops=dict(
            arrowstyle="-|>",
            lw=1.2,
            color=GUIDE,
            shrinkA=5,
            shrinkB=4,
            connectionstyle=f"arc3,rad={rad}",
        ),
    )


plt.rcParams.update(
    {
        "figure.facecolor": FIG_BG,
        "axes.facecolor": FIG_BG,
        "savefig.facecolor": FIG_BG,
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "axes.titleweight": 600,
        "axes.titlesize": 24,
    }
)

fig = plt.figure(figsize=(15.5, 8.45))
gs = fig.add_gridspec(3, 1, height_ratios=[0.85, 3.0, 2.55], hspace=0.18)
header_ax = fig.add_subplot(gs[0])
ax = fig.add_subplot(gs[1])
legend_ax = fig.add_subplot(gs[2])

header_ax.axis("off")
header_ax.set_xlim(0, 1)
header_ax.set_ylim(0, 1)
header_ax.text(
    0.5,
    0.78,
    "SAR Lamination Ladder Mixer",
    ha="center",
    va="center",
    fontsize=30,
    fontweight=600,
    color=TEXT,
    transform=header_ax.transAxes,
)
header_ax.text(
    0.5,
    0.34,
    "2D unit-cell sketch with vertically stretched channel for presentation clarity",
    ha="center",
    va="center",
    fontsize=15,
    color="#48677e",
    transform=header_ax.transAxes,
)

ax.set_aspect("equal", adjustable="box")

channel = Rectangle((0, 0), L, H_plot, facecolor="white", edgecolor=CHANNEL_EDGE, linewidth=3.2)
ax.add_patch(channel)

for k in range(N):
    xk = L0 + k * L_cell

    ax.plot([xk, xk], [0, H_plot], color=GUIDE_SOFT, lw=1.2, zorder=0)
    ax.plot([xk + L_cell, xk + L_cell], [0, H_plot], color=GUIDE_SOFT, lw=1.2, zorder=0)

    xs0, xs1 = xk, xk + L_s
    y0s, y1s = (H - t_s) / 2, (H + t_s) / 2
    ax.add_patch(
        Rectangle(
            (xs0, ys(y0s)),
            L_s,
            ys(t_s),
            facecolor=SPLIT_FACE,
            edgecolor=SPLIT_EDGE,
            linewidth=2.0,
        )
    )

    xc0, xc1 = xk + L_s, xk + L_s + L_c
    xx = np.linspace(xc0, xc1, 400)
    xi = xx - xc0
    bump = h_d * g_env(xi, L_c)
    top_curve = H - np.clip(bump + delta, 0, H)

    ax.fill_between(xx, 0, ys(bump), color=BOTTOM_FACE, alpha=0.95, linewidth=0)
    ax.plot(xx, ys(bump), color=BOTTOM_EDGE, linewidth=2.1)

    ax.fill_between(xx, ys(top_curve), H_plot, color=TOP_FACE, alpha=0.92, linewidth=0)
    ax.plot(xx, ys(top_curve), color=TOP_EDGE, linewidth=2.1)

    xm0, xm1 = xk + L_s + L_c, xk + L_cell
    y0m, y1m = (H - t_m) / 2, (H + t_m) / 2
    ax.add_patch(
        Rectangle(
            (xm0, ys(y0m)),
            L_m,
            ys(t_m),
            facecolor=MERGE_FACE,
            edgecolor=MERGE_EDGE,
            linewidth=2.0,
        )
    )

repeat_y = H_plot + 1.12
ax.plot([L0, L0 + N * L_cell], [repeat_y, repeat_y], color=GUIDE, lw=1.4)
ax.plot([L0, L0], [repeat_y - 0.09, repeat_y + 0.09], color=GUIDE, lw=1.4)
ax.plot([L0 + N * L_cell, L0 + N * L_cell], [repeat_y - 0.09, repeat_y + 0.09], color=GUIDE, lw=1.4)
ax.text(
    L0 + 0.5 * N * L_cell,
    repeat_y + 0.18,
    rf"$N = {N}$ repeated unit cells",
    ha="center",
    va="bottom",
    fontsize=12,
    color=TEXT,
    bbox=dict(boxstyle="round,pad=0.22", fc=PANEL_BG, ec="none"),
)

ax.plot([0, 0], [H_plot, H_plot + 0.42], color=GUIDE, lw=1.1)
ax.plot([L0, L0], [H_plot, H_plot + 0.42], color=GUIDE, lw=1.1)
dim_h(ax, 0, L0, H_plot + 0.42, r"$L_0$")

ax.plot([L0, L0], [0, -0.34], color=GUIDE, lw=1.1)
ax.plot([L0 + L_cell, L0 + L_cell], [0, -0.34], color=GUIDE, lw=1.1)
dim_h(ax, L0, L0 + L_cell, -0.34, r"$L_{\mathrm{cell}}$", text_dy=0.06)

ax.plot([0, 0], [0, -1.10], color=GUIDE, lw=1.1)
ax.plot([L, L], [0, -1.10], color=GUIDE, lw=1.1)
dim_h(ax, 0, L, -1.10, r"$L$", text_dy=0.09)

ax.plot([-0.86, 0], [0, 0], color=GUIDE, lw=1.1)
ax.plot([-0.86, 0], [H_plot, H_plot], color=GUIDE, lw=1.1)
dim_v(ax, -0.86, 0, H_plot, r"$H$")

callout(ax, r"$\Gamma_{\mathrm{in}}$", (0, 0.5 * H_plot), (-1.45, 0.54 * H_plot), ha="right", rad=0.08)
callout(ax, r"$\Gamma_{\mathrm{out}}$", (L, 0.5 * H_plot), (L + 1.45, 0.54 * H_plot), ha="left", rad=-0.08)

first_cell_x = L0
split_mid_x = first_cell_x + 0.45 * L_s
interaction_mid_x = first_cell_x + L_s + 0.5 * L_c
merge_mid_x = first_cell_x + L_s + L_c + 0.55 * L_m

callout(
    ax,
    r"$t_s$",
    (split_mid_x, ys((H + t_s) / 2)),
    (first_cell_x + 0.45, H_plot + 0.78),
    rad=-0.12,
)
callout(
    ax,
    r"$w_s$",
    (split_mid_x, ys(0.25 * ((H - t_s) / 2))),
    (first_cell_x + 0.90, -0.48),
    rad=0.12,
)
callout(
    ax,
    r"$h_d$",
    (interaction_mid_x, ys(h_d)),
    (first_cell_x + L_s + 0.18 * L_c, -0.88),
    rad=-0.10,
)
callout(
    ax,
    r"$\delta$",
    (interaction_mid_x, ys(H - h_d - delta)),
    (first_cell_x + L_s + 0.60 * L_c, H_plot + 0.98),
    rad=0.12,
)
callout(
    ax,
    r"$t_m$",
    (merge_mid_x, ys((H + t_m) / 2)),
    (first_cell_x + L_s + L_c + 0.18 * L_m, H_plot + 0.74),
    rad=-0.08,
)

ax.set_xlim(-1.65, L + 1.65)
ax.set_ylim(-1.34, H_plot + 1.52)
ax.axis("off")

legend_ax.axis("off")
legend_ax.set_xlim(0, 1)
legend_ax.set_ylim(0, 1)
legend_ax.add_patch(
    FancyBboxPatch(
        (0.012, 0.05),
        0.976,
        0.90,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        facecolor=PANEL_BG,
        edgecolor=GUIDE_SOFT,
        linewidth=1.2,
        transform=legend_ax.transAxes,
    )
)

legend_handles = [
    Patch(facecolor=SPLIT_FACE, edgecolor=SPLIT_EDGE, label="split splitter"),
    Patch(facecolor=BOTTOM_FACE, edgecolor=BOTTOM_EDGE, label="bottom deflector"),
    Patch(facecolor=TOP_FACE, edgecolor=TOP_EDGE, label="top deflector"),
    Patch(facecolor=MERGE_FACE, edgecolor=MERGE_EDGE, label="merge splitter"),
]
legend_ax.legend(
    handles=legend_handles,
    loc="upper left",
    bbox_to_anchor=(0.03, 0.94),
    ncol=4,
    frameon=False,
    fontsize=11,
    handlelength=1.8,
    columnspacing=1.8,
)

legend_ax.text(0.03, 0.70, "Symbol key", fontsize=13, fontweight="bold", color=TEXT, transform=legend_ax.transAxes)
legend_ax.text(0.53, 0.70, "Parameter key", fontsize=13, fontweight="bold", color=TEXT, transform=legend_ax.transAxes)

symbol_text = (
    r"$\Gamma_{\mathrm{in}}$: inlet boundary" "\n"
    r"$\Gamma_{\mathrm{out}}$: outlet boundary" "\n"
    r"$H$: channel height" "\n"
    r"$L$: total mixer length" "\n"
    r"$L_0$: inlet/outlet buffer length" "\n"
    r"$L_{\mathrm{cell}}$: SAR unit-cell length" "\n"
    rf"$N$: number of repeated unit cells ({N})"
)
param_text = (
    r"$t_s$: split-section splitter thickness" "\n"
    r"$t_m$: merge-section splitter thickness" "\n"
    r"$w_s$: nominal subchannel width after split" "\n"
    r"$h_d$: deflector intrusion height from wall" "\n"
    r"$\delta$: top-deflector vertical bias" "\n"
    r"$L_s, L_c, L_m$: split, interaction, merge lengths" "\n"
    r"$r$: fillet radius design parameter"
)

legend_ax.text(0.03, 0.64, symbol_text, fontsize=10.4, color=TEXT, va="top", transform=legend_ax.transAxes)
legend_ax.text(0.53, 0.64, param_text, fontsize=10.4, color=TEXT, va="top", transform=legend_ax.transAxes)

fig.text(
    0.5,
    0.028,
    rf"Display note: channel height shown with a {Y_SCALE:.2f}x y-scale stretch for readability.",
    ha="center",
    va="center",
    fontsize=10.2,
    color="#60798c",
    bbox=dict(boxstyle="round,pad=0.28", fc="#f0f5fa", ec="none"),
)

png_path = Path(__file__).with_suffix(".png")
fig.savefig(png_path, dpi=300, bbox_inches="tight")
print(f"Saved {png_path}")
