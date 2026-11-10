#!/usr/bin/env python3
"""Band-structure figure from the bands.x output (runs/01_scf_bands/BN.bands.dat.gnu).

Energies are shifted so that the valence-band maximum (VBM) is 0 eV; the k-axis is the
path length in 1/Angstrom. Occupied bands (1-4) solid blue, unoccupied (5-8) dashed orange.
Writes runs/bands.pdf (for LaTeX) and runs/bands.png.
"""
import pathlib
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RUNS = pathlib.Path(__file__).resolve().parents[1]
D = RUNS / "01_scf_bands"
A = 2.504                                   # Angstrom (BN.in)
TPIBA = 2 * np.pi / A                       # bands.x x-axis is in units of 2pi/a
NOCC = 4                                    # 8 valence electrons / 2

gnu = np.loadtxt(D / "BN.bands.dat.gnu")
nk = int(np.argmax(np.diff(gnu[:, 0]) < 0) + 1) if np.any(np.diff(gnu[:, 0]) < 0) else len(gnu)
x = gnu[:nk, 0] * TPIBA
E = gnu[:, 1].reshape(-1, nk)               # (nbnd, nk)
vbm = E[NOCC - 1].max()
E = E - vbm
cbm = E[NOCC].min()
iv, ic = np.argmax(E[NOCC - 1]), np.argmin(E[NOCC])

# high-symmetry positions printed by bands.x
hs = [float(m) for m in re.findall(r"high-symmetry point:.*x coordinate\s+([-\d.]+)",
                                    (D / "bands_pp.out").read_text())]
hs = np.array(hs) * TPIBA
labels = [r"$\Gamma$", r"$K$", r"$M$"]

BLUE, ORANGE = "#2a78d6", "#eb6834"         # validated categorical slots 1-2 (light surface)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"

plt.rcParams.update({"font.size": 9, "font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(figsize=(4.3, 3.3))
ax.axhspan(0, cbm, color="#f0efec", zorder=0, lw=0)          # band gap region
for xv in hs[1:-1]:
    ax.axvline(xv, color=GRID, lw=0.8, zorder=1)
ax.axhline(0, color=MUTED, lw=0.6, ls=":", zorder=1)

for n in range(E.shape[0]):
    occ = n < NOCC
    ax.plot(x, E[n], color=BLUE if occ else ORANGE, ls="-" if occ else (0, (4, 2)),
            lw=1.4, zorder=3, label=(None if n not in (0, NOCC) else
                                     ("occupied (bands 1–4)" if occ else "unoccupied (bands 5–8)")))

# band edges + gap annotation
ax.plot([x[iv]], [0], "o", ms=5, color=BLUE, mec="white", mew=1.0, zorder=5)
ax.plot([x[ic]], [cbm], "o", ms=5, color=ORANGE, mec="white", mew=1.0, zorder=5)
ax.annotate("", xy=(x[ic], cbm), xytext=(x[iv], 0),
            arrowprops=dict(arrowstyle="<->", color=INK, lw=0.8, shrinkA=3, shrinkB=3), zorder=6)
ax.text(x[iv] + 0.05, cbm / 2, f"$E_g$ = {cbm:.2f} eV\n(direct, at $K$)",
        va="center", ha="left", fontsize=8, color=INK, zorder=6,
        bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.85))

ax.set_xlim(hs[0], hs[-1])
ax.set_ylim(-19, 12)
ax.set_xticks(hs, labels)
ax.set_ylabel(r"$\varepsilon_{n\mathbf{k}} - \varepsilon_{\mathrm{VBM}}$ (eV)")
ax.set_xlabel(r"wave vector $\mathbf{k}$ (path drawn to scale in Å$^{-1}$)", fontsize=8, color=MUTED)
ax.tick_params(axis="x", length=0, labelsize=10, labelcolor=INK)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="lower right", fontsize=7.5, frameon=True, framealpha=1, edgecolor=GRID)
fig.tight_layout(pad=0.3)
fig.savefig(RUNS / "bands.pdf")
fig.savefig(RUNS / "bands.png", dpi=200)
print(f"VBM index {iv} (x={x[iv]:.4f}), CBM index {ic} (x={x[ic]:.4f}), gap {cbm:.4f} eV, "
      f"high-symmetry x = {np.round(hs, 4)} 1/A")
