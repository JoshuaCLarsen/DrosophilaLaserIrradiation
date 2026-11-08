#!/usr/bin/env python3
"""
Monte Carlo dose model runs for article_revised_monte_6 (physics in freefall_dose_model_monte.py).

Two position models, each a population of N_fly simulated flies that fall straight down through 150 planes
(1 ms steps, t_exp = 150 ms) at a fixed transverse position:
  uniform : entry point uniform over the accessible disc (baseline)
  band    : entry point uniform over the band within one body radius of the stem wall
Only the dose-response parameters k and F50 (Eq. 8) are fitted; every curve approaches unity at high power.

Outputs:  lethalityErrorBars.txt, monte6Params.txt, monte6Curves.txt, monte6Energy.npz in this directory
"""
import numpy as np
import os
from freefall_dose_model_monte import (simulate, fit, mortality, logistic, N_exp, n_per, A_fly, seed, OUTPUT_DIR, HERE,
                                       z_top, dz, N_z, R_acc, r_f)

data = np.loadtxt(os.path.join(HERE, "flyDeathRates.txt"), delimiter=",")
powers, deaths = data[:, 0], data[:, 1]
Pmax = powers.max()

# %% +/- 1 sigma binomial error bars of the measured lethality (read by the pgfplots figure)
with open(os.path.join(OUTPUT_DIR, "lethalityErrorBars.txt"), "w") as fe:
    fe.write("power,measured,se\n")
    for P, D in zip(powers, deaths):
        fe.write(f"{P:.4f},{D:.5f},{np.sqrt(D * (1 - D) / n_per):.4f}\n")

# %% simulate and fit
res = {}
for name in ("uniform", "band"):
    E1, track = simulate(name, 0.0, 0.0, 0.0, n_track=60)
    f = fit(E1, powers, deaths)
    res[name] = dict(E1=E1, track=track, f=f)
    F = E1 * Pmax / (A_fly / 100)
    print(f"[{name}] E1/W mean {E1.mean():.4f} J, p10-p90 {np.percentile(E1,10):.4f}-{np.percentile(E1,90):.4f}, max {E1.max():.4f}, "
          f"CV {E1.std()/E1.mean():.2f}; mean F at Pmax {F.mean():.2f} J/cm2")
    print(f"   k={f['k']:.3f} F50={f['F50']:.3f} LD90={f['LD90']:.2f} E90={f['E90']:.3f} R2={f['R2']:.3f} RMSE={f['RMSE']:.3f} "
          f"chi2={f['chisq']:.1f}/{f['dof']} p={f['p']:.3f} M(Pmax)={f['pred'][-1]:.3f} below F50 at Pmax {np.mean(F < f['F50'])*100:.1f}%")

with open(os.path.join(OUTPUT_DIR, "monte6Params.txt"), "w") as fh:
    fh.write("model,k,F50,LD90,E90,R2,RMSE,chisq,dof,p,meanE1,cvE1\n")
    for name, r in res.items():
        f = r["f"]
        fh.write(f"{name},{f['k']:.4f},{f['F50']:.4f},{f['LD90']:.4f},{f['E90']:.4f},{f['R2']:.4f},{f['RMSE']:.4f},"
                 f"{f['chisq']:.3f},{f['dof']},{f['p']:.4f},{r['E1'].mean():.5f},{r['E1'].std()/r['E1'].mean():.4f}\n")

# %% curves and the +/- 1 sigma scatter of simulated 100-fly experiments for the baseline
rng = np.random.default_rng(seed + 21)
Pfine = np.linspace(0, 12, 121)
E1b, fb = res["uniform"]["E1"], res["uniform"]["f"]
with open(os.path.join(OUTPUT_DIR, "monte6Curves.txt"), "w") as fc:
    fc.write("power,uniform,uniform_lo,uniform_hi,band\n")
    for P in Pfine:
        pick = rng.integers(0, len(E1b), (N_exp, n_per))
        dead = rng.uniform(size=(N_exp, n_per)) < logistic(E1b[pick] * P / (A_fly / 100), fb["k"], fb["F50"])
        sd = dead.mean(axis=1).std()
        m = {k_: mortality(r["E1"], P, r["f"]["k"], r["f"]["F50"]) for k_, r in res.items()}
        fc.write(f"{P:.3f},{m['uniform']:.5f},{m['uniform']-sd:.5f},{m['uniform']+sd:.5f},{m['band']:.5f}\n")
        if abs(P - 12) < 1e-9 or abs(P - 7.1) < 1e-9:
            print(f"   P={P:.1f} W: uniform {m['uniform']:.3f}, band {m['band']:.3f}; SD {sd:.3f} vs binomial {np.sqrt(m['uniform']*(1-m['uniform'])/n_per):.3f}")

np.savez(os.path.join(HERE, "monte6Energy.npz"), **{f"E1_{k_}": r["E1"] for k_, r in res.items()},
         **{f"track_{k_}": r["track"] for k_, r in res.items()},
         **{f"LD90_{k_}": r["f"]["LD90"] for k_, r in res.items()}, z_top=z_top, dz=dz, N_z=N_z, R_acc=R_acc, r_f=r_f)
print("wrote monte6Params.txt, monte6Curves.txt, monte6Energy.npz")
