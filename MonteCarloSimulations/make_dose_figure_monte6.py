#!/usr/bin/env python3
"""Dose figure for article_revised_monte_6: monte6Dose.pdf in this directory
(a) radial position of the fly centre versus depth for sample Monte Carlo flies of the two position models,
    with the 1/e half-widths of the beam; (b) distribution of the body-averaged fluence per fly at the maximum power.
Run freefall_dose_model_monte6.py first."""
import numpy as np, os, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from freefall_dose_model_monte import w0x, w0y, thx, thy, R_tube, stem_len, A_fly, OUTPUT_DIR, HERE

d = np.load(os.path.join(HERE, "monte6Energy.npz"))
Pmax = np.loadtxt(os.path.join(HERE, "flyDeathRates.txt"), delimiter=",")[:, 0].max()
depth = (np.arange(int(d["N_z"])) + 0.5) * float(d["dz"])
z = float(d["z_top"]) + depth
names = [("uniform", "Uniform", "#2a78d6", "-"),
         ("band", "Wall band", "#eb6834", "--")]

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "axes.edgecolor": "#6b6b6b", "xtick.color": "#4a4a4a",
                     "ytick.color": "#4a4a4a", "axes.labelcolor": "#262626"})
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.6), gridspec_kw={"width_ratios": [1, 1.3]})

# (a) radial position versus depth
a = ax[0]
a.fill_betweenx(depth, 0, np.minimum(w0x + thx * z, R_tube), color="#c6dbef", lw=0, label="Beam $1/e$ half-width ($x$ to $y$)")
a.fill_betweenx(depth, np.minimum(w0x + thx * z, R_tube), np.minimum(w0y + thy * z, R_tube), color="#deebf7", lw=0)
a.plot([R_tube, R_tube], [0, stem_len], color="#4a4a4a", lw=1.0)
a.plot([float(d["R_acc"])] * 2, [0, stem_len], color="#4a4a4a", lw=0.5, ls=":")
rng = np.random.default_rng(5)
for key, lab, c, ls in names:
    tr = d[f"track_{key}"]
    r = np.hypot(tr[:, 0, :], tr[:, 1, :])
    pick = rng.choice(r.shape[1], 5, replace=False)
    for i, j in enumerate(pick):
        a.plot(r[:, j], depth, color=c, ls=ls, lw=1.1, label=lab if i == 0 else None)
a.axhline(stem_len, color="#9e9e9e", lw=0.5, ls=":")
a.text(0.05, stem_len - 1.0, "stem exit", fontsize=7, color="#6b6b6b")
a.set_ylim(depth[-1] + 0.5, 0); a.set_xlim(0, 3.3)
a.set_xlabel("Radial position of fly center (mm)"); a.set_ylabel("Depth below stem entrance (mm)")
a.set_title("(a) Simulated fly paths", fontsize=9, loc="left")

# (b) fluence distribution at the maximum power
b = ax[1]
edges = np.linspace(0, 13.5, 55)
hs = {key: np.histogram(d[f"E1_{key}"] * Pmax / (A_fly / 100), bins=edges)[0] / len(d[f"E1_{key}"]) for key, *_ in names}
top = max(h.max() for h in hs.values())
for key, lab, c, ls in names:
    b.step(edges[:-1], hs[key], where="post", color=c, ls=ls, lw=1.6, label=lab)
    b.plot([float(d[f"LD90_{key}"])] * 2, [1.08 * top, 1.16 * top], color=c, lw=2.2, solid_capstyle="round")
b.text(1.5, 1.12 * top, "$F_{90}$", fontsize=8, color="#4a4a4a", va="center")
b.set_xlim(0, 13.5); b.set_ylim(0, 1.2 * top)
b.set_xlabel(f"Body-averaged fluence at {Pmax:.2f} W (J/cm$^2$)"); b.set_ylabel("Fraction of flies")
b.set_title(f"(b) Dose delivered at {Pmax:.2f} W", fontsize=9, loc="left")
b.grid(axis="y", color="#e5e5e5", lw=0.5)
for a_ in ax: a_.spines[["top", "right"]].set_visible(False)
h1, l1 = ax[0].get_legend_handles_labels()
fig.legend(h1, l1, loc="lower center", ncol=3, frameon=False, fontsize=7.5, handlelength=3)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig(os.path.join(OUTPUT_DIR, "monte6Dose.pdf")); fig.savefig(os.path.join(OUTPUT_DIR, "monte6Dose.png"), dpi=150)
print("wrote monte6Dose.pdf")
