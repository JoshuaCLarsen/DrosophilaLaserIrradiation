#!/usr/bin/env python3
"""
Monte Carlo free-fall dose model for Drosophila hydei irradiated inside a funnel stem.

Each simulated fly is a disc of area A_f that falls through the 50 mm interaction length (40 mm stem + 10 mm gap
to the slide) in N_z = 150 planes of thickness dz = L/N_z, i.e. time steps of dt = t_exp/N_z = 1 ms.  In every
plane the fly collects the beam power that falls on its disc for dt, and then moves to the next plane:

    x_{i+1} = x_i + vx dt,    y_{i+1} = y_i + vy dt,    z_{i+1} = z_i + dz

so that its transverse position is correlated from plane to plane through its lateral velocity (vx, vy).
A fly reaching the stem wall loses the normal component of its lateral velocity (coefficient of restitution e;
e = 0: the fly stops at the wall and slides down it, e = 1: specular rebound).

Entry condition ("cone"): the fly slides down the 60-degree funnel cone and arrives at the stem wall at a random
azimuth, moving inward along the fall line of the cone.  Its lateral speed follows from the cone half-angle alpha
and the mean vertical speed v_z = L / t_exp:  v_lat = v_z tan(alpha).
Entry condition ("band"): the fly centre is uniformly distributed over the band R_acc - r_f <= r <= R_acc next to
the wall (the "wall" position distribution of the earlier dose model).
Entry condition ("uniform", v_lat = 0): the fly enters at a uniformly random position and falls straight down;
this reproduces the fixed-position integral of the earlier dose model (Model/freefall_dose_model.py).

The energy collected by fly j at 1 W is E1_j (J/W).  Because the irradiance is proportional to the laser power, the
energy at power P is E1_j P and the body-averaged fluence is F_j = E1_j P / A_f.  Each fly dies with probability
p(F_j) = 1 / (1 + exp(-k (F_j - F50)))  (Eq. 8), and the predicted mortality is the fraction of dead flies,
M(P) = (1/N) sum_j p(F_j).  k and F50 are fitted to the ten measured mortality fractions.

Outputs (written to this script's directory, suffix "Monte"):
  monteModelParams.txt     fitted parameters and goodness of fit for each configuration
  monteModelCurves.txt     predicted mortality versus power (+/- 1 sigma of simulated 100-fly experiments)
  monteDosePercentiles.txt per-watt energy statistics of the simulated population
  monteEnergy.npz          per-fly energies and sample trajectories (read by make_dose_figure_monte.py)

Usage:  python3 freefall_dose_model_monte.py            (baseline + comparison configurations)
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gamma
from scipy.stats import chi2
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = HERE

# %% ------------------------------------------------------------------ experimental parameters (Table 3)
f_lens    = 45.0                 # mm, focal length of the module lens
gap_laser = 31.0                 # mm, lens to top of funnel cone
cone_len  = 30.0                 # mm, cone height
stem_len  = 40.0                 # mm, stem length
slide_gap = 10.0                 # mm, stem exit to microscope slide
R_tube    = 3.0                  # mm, stem inner radius (6 mm ID)
alpha     = np.deg2rad(30.0)     # cone half-angle (60-degree analytical funnel)
t_exp     = 0.150                # s, measured transit time over the interaction length
A_fly     = 3.3                  # mm^2, projected area of an adult D. hydei
# beam: waist = 0.08 x 0.10 mm datasheet focal spot; 1/e half-widths at 1.0 m from the object-space fit (Table 1)
w0x, w0y  = 0.040, 0.050         # mm
W1x, W1y  = 36.81, 40.58         # mm, at z1 = 955 mm beyond the focal plane
z1        = 1000.0 - f_lens
nx, ny    = 2, 6                 # super-Gaussian exponents

# %% ------------------------------------------------------------------ simulation parameters
N_fly  = 20000                   # simulated flies per configuration
N_z    = 150                     # planes across the interaction length (dt = 1 ms)
N_disc = 61                      # sample points on the fly disc (Vogel spiral, equal area per point)
N_exp  = 2000                    # simulated experiments of n = 100 flies (error band)
n_per  = 100                     # flies per power setting in the experiment
seed   = 1

L      = stem_len + slide_gap
z_top  = gap_laser + cone_len - f_lens          # 16 mm beyond the focal plane
dz     = L / N_z
dt     = t_exp / N_z
v_z    = L / t_exp                               # mm/s
r_f    = np.sqrt(A_fly / np.pi)                  # mm, radius of the fly disc
R_acc  = R_tube - r_f                            # mm, largest radius of the fly centre inside the stem
v_lat_cone = v_z * np.tan(alpha)                 # mm/s

# %% ------------------------------------------------------------------ beam (Eqs. 2, 6)
thx = (W1x - w0x) / z1
thy = (W1y - w0y) / z1
Nrm = 4 * gamma(1 + 1 / nx) * gamma(1 + 1 / ny)

def irradiance(x, y, z):
    """Irradiance (W/mm^2) per watt of laser power at transverse position (x, y), distance z beyond the focus."""
    wx, wy = w0x + thx * z, w0y + thy * z
    return np.exp(-np.abs(x / wx) ** nx - np.abs(y / wy) ** ny) / (Nrm * wx * wy)

# fly footprint: N_disc points spread evenly over a disc of radius r_f (each point represents A_fly / N_disc)
kk = np.arange(N_disc) + 0.5
disc_x = r_f * np.sqrt(kk / N_disc) * np.cos(kk * np.pi * (3 - np.sqrt(5)))
disc_y = r_f * np.sqrt(kk / N_disc) * np.sin(kk * np.pi * (3 - np.sqrt(5)))

# %% ------------------------------------------------------------------ Monte Carlo fall
def simulate(entry="cone", v_lat=v_lat_cone, e=0.0, spread_deg=0.0, n=N_fly, n_track=40, rng=None):
    """Drop n flies through the interaction length.  Returns E1 (J per W of laser power, one value per fly) and the
    transverse positions of the first n_track flies at every plane (for plotting)."""
    rng = np.random.default_rng(seed) if rng is None else rng
    phi = rng.uniform(0, 2 * np.pi, n)
    if entry == "cone":                          # at the wall, moving inward along the fall line of the cone
        x, y = R_acc * np.cos(phi), R_acc * np.sin(phi)
        heading = phi + np.pi + np.deg2rad(spread_deg) * rng.uniform(-1, 1, n)
    elif entry == "band":                        # uniformly distributed over the band within one body radius of the wall
        r = np.sqrt(rng.uniform((R_acc - r_f) ** 2, R_acc ** 2, n))
        x, y = r * np.cos(phi), r * np.sin(phi)
        heading = rng.uniform(0, 2 * np.pi, n)
    else:                                        # uniformly distributed over the accessible disc
        r = R_acc * np.sqrt(rng.uniform(0, 1, n))
        x, y = r * np.cos(phi), r * np.sin(phi)
        heading = rng.uniform(0, 2 * np.pi, n)
    vx, vy = v_lat * np.cos(heading), v_lat * np.sin(heading)

    E1 = np.zeros(n)
    track = np.zeros((N_z, 2, n_track))
    for i in range(N_z):
        z = z_top + (i + 0.5) * dz
        track[i] = x[:n_track], y[:n_track]
        # energy collected in this plane: mean irradiance over the disc x fly area x time step
        I = irradiance(x[:, None] + disc_x, y[:, None] + disc_y, z)
        E1 += I.mean(axis=1) * A_fly * dt
        # move to the next plane
        x = x + vx * dt
        y = y + vy * dt
        # stem wall: keep the fly centre inside R_acc and remove (e = 0) or reverse (e = 1) its normal velocity
        if z - z_top < stem_len:
            r = np.hypot(x, y)
            hit = r > R_acc
            if hit.any():
                ux, uy = x[hit] / r[hit], y[hit] / r[hit]
                x[hit], y[hit] = R_acc * ux, R_acc * uy
                vn = vx[hit] * ux + vy[hit] * uy
                vx[hit] -= (1 + e) * vn * ux
                vy[hit] -= (1 + e) * vn * uy
    return E1, track

# %% ------------------------------------------------------------------ mortality and fitting (Eq. 8)
def logistic(F, k, F50):
    return 1.0 / (1.0 + np.exp(-k * (F - F50)))

def mortality(E1, P, k, F50, s=0.0):
    """Fraction of the simulated population that dies at laser power P."""
    F = E1 * P / (A_fly / 100.0)                 # J/cm^2
    return (1 - s) * logistic(F, k, F50).mean()

def fit(E1, powers, deaths, with_s=False):
    def resid(q):
        s = q[2] if with_s else 0.0
        return np.array([mortality(E1, P, q[0], q[1], s) for P in powers]) - deaths
    q0, lb, ub = ([1.0, 4.0, 0.1], [1e-3, 0, 0], [20, 50, 1]) if with_s else ([1.0, 4.0], [1e-3, 0], [20, 50])
    q = least_squares(resid, q0, bounds=(lb, ub)).x
    k, F50, s = q[0], q[1], (q[2] if with_s else 0.0)
    pred = np.array([mortality(E1, P, k, F50, s) for P in powers])
    ss_res = np.sum((deaths - pred) ** 2)
    # chi-square against the binomial variance of the measured fractions (control point excluded: zero variance)
    use = (deaths > 0) & (deaths < 1)
    chisq = np.sum((deaths[use] - pred[use]) ** 2 / (deaths[use] * (1 - deaths[use]) / n_per))
    dof = use.sum() - len(q)
    # log-likelihood of the observed counts (for likelihood-ratio tests)
    pc = np.clip(pred, 1e-9, 1 - 1e-9)
    lnL = np.sum(n_per * (deaths * np.log(pc) + (1 - deaths) * np.log(1 - pc)))
    LD90 = F50 + np.log(9) / k
    return dict(k=k, F50=F50, s=s, pred=pred, LD90=LD90, E90=LD90 * A_fly / 100, R2=1 - ss_res / np.sum((deaths - deaths.mean()) ** 2),
                RMSE=np.sqrt(ss_res / len(deaths)), chisq=chisq, dof=dof, p=chi2.sf(chisq, dof), lnL=lnL)

# %% ------------------------------------------------------------------ run
if __name__ == "__main__":
    data = np.loadtxt(os.path.join(HERE, "flyDeathRates.txt"), delimiter=",")
    powers, deaths = data[:, 0], data[:, 1]
    Pmax = powers.max()
    print(f"v_z = {v_z:.0f} mm/s, cone half-angle {np.rad2deg(alpha):.0f} deg -> v_lat = {v_lat_cone:.0f} mm/s; "
          f"r_f = {r_f:.2f} mm, accessible radius {R_acc:.2f} mm; dz = {dz:.3f} mm, dt = {dt*1e3:.1f} ms")
    print(f"a fly entering at the wall crosses the accessible disc in {2*R_acc/v_lat_cone*1e3:.0f} ms ({2*R_acc/v_lat_cone*v_z:.1f} mm of fall)")

    configs = {                                   # name: (entry, v_lat, e, heading spread in degrees)
        "cone":       ("cone", v_lat_cone, 0.0, 0.0),    # baseline: slide off the cone, cross the bore, slide down the far wall
        "coneBounce": ("cone", v_lat_cone, 1.0, 0.0),    # sensitivity: specular rebound at the wall
        "coneSpread": ("cone", v_lat_cone, 0.0, 90.0),   # sensitivity: any inward heading at entry
        "straight":   ("uniform", 0.0, 0.0, 0.0),        # fixed transverse position (earlier integral model)
    }
    out = {}
    fh = open(os.path.join(OUTPUT_DIR, "monteModelParams.txt"), "w")
    fh.write("config,fit,k,F50,s,LD90,E90,R2,RMSE,chisq,dof,p,LRp_s,meanE1,cvE1,fracBelowF50atPmax,M_at_Pmax\n")
    for name, (entry, vl, e, sp) in configs.items():
        E1, track = simulate(entry, vl, e, sp)
        f2 = fit(E1, powers, deaths)
        f3 = fit(E1, powers, deaths, with_s=True)
        LRp = chi2.sf(max(2 * (f3["lnL"] - f2["lnL"]), 0), 1)
        F1 = E1 / (A_fly / 100)
        out[name] = dict(E1=E1, track=track, f2=f2, f3=f3)
        print(f"\n[{name}] E1 per watt: mean {E1.mean():.4f} J, p10 {np.percentile(E1,10):.4f}, median {np.median(E1):.4f}, "
              f"p90 {np.percentile(E1,90):.4f}, max {E1.max():.4f}  (CV {E1.std()/E1.mean():.2f})")
        for lab, f in (("2", f2), ("3", f3)):
            below = np.mean(F1 * Pmax < f["F50"])
            print(f"   fit {lab} param: k={f['k']:.3f} F50={f['F50']:.3f} s={f['s']:.3f} LD90={f['LD90']:.2f} J/cm2 E90={f['E90']:.3f} J "
                  f"R2={f['R2']:.3f} RMSE={f['RMSE']:.3f} chi2={f['chisq']:.1f}/{f['dof']} p={f['p']:.3f} M(Pmax)={f['pred'][-1]:.3f}; "
                  f"{below*100:.1f}% of flies below F50 at {Pmax} W")
            fh.write(f"{name},{lab},{f['k']:.4f},{f['F50']:.4f},{f['s']:.4f},{f['LD90']:.4f},{f['E90']:.4f},{f['R2']:.4f},{f['RMSE']:.4f},"
                     f"{f['chisq']:.3f},{f['dof']},{f['p']:.4f},{LRp:.2e},{E1.mean():.5f},{E1.std()/E1.mean():.4f},{below:.4f},{f['pred'][-1]:.4f}\n")
        print(f"   likelihood-ratio test for s: p = {LRp:.2e}")
    fh.close()

    # %% simulated experiments: n = 100 flies per power setting, each fly a non-responder with probability s
    base = out["cone"]; E1 = base["E1"]; f3 = base["f3"]; f2 = base["f2"]
    rng = np.random.default_rng(seed + 1)
    Pfine = np.linspace(0, 8, 81)
    with open(os.path.join(OUTPUT_DIR, "monteModelCurves.txt"), "w") as fc:
        fc.write("power,coneS,coneS_lo,coneS_hi,cone2,straight\n")
        for P in Pfine:
            pick = rng.integers(0, len(E1), (N_exp, n_per))                     # flies drawn from the simulated population
            responder = rng.uniform(size=(N_exp, n_per)) >= f3["s"]
            dead = responder & (rng.uniform(size=(N_exp, n_per)) < logistic(E1[pick] * P / (A_fly / 100), f3["k"], f3["F50"]))
            frac = dead.mean(axis=1)
            m = mortality(E1, P, f3["k"], f3["F50"], f3["s"])
            fc.write(f"{P:.3f},{m:.5f},{m-frac.std():.5f},{m+frac.std():.5f},{mortality(E1, P, f2['k'], f2['F50']):.5f},"
                     f"{mortality(out['straight']['E1'], P, out['straight']['f2']['k'], out['straight']['f2']['F50']):.5f}\n")
            if np.any(np.isclose(P, powers, atol=0.05)):
                print(f"   P={P:.1f} W: model {m:.3f}, SD of simulated 100-fly experiments {frac.std():.3f}, binomial {np.sqrt(m*(1-m)/n_per):.3f}")

    # %% per-watt dose statistics and data for the figures
    with open(os.path.join(OUTPUT_DIR, "monteDosePercentiles.txt"), "w") as fd:
        fd.write("config,mean,p10,p50,p90,min,max\n")
        for name in out:
            e1 = out[name]["E1"]
            fd.write(f"{name},{e1.mean():.5f},{np.percentile(e1,10):.5f},{np.median(e1):.5f},{np.percentile(e1,90):.5f},{e1.min():.5f},{e1.max():.5f}\n")
    np.savez(os.path.join(HERE, "monteEnergy.npz"), E1_cone=E1, E1_straight=out["straight"]["E1"], track=base["track"],
             k=f3["k"], F50=f3["F50"], s=f3["s"], LD90=f3["LD90"], z_top=z_top, dz=dz, N_z=N_z, R_acc=R_acc, r_f=r_f)
    print(f"\nbaseline: fly centre crosses the accessible disc in {2*R_acc/v_lat_cone*1e3:.0f} ms; "
          f"mean body-averaged fluence at {Pmax} W = {E1.mean()*Pmax/(A_fly/100):.2f} J/cm2 "
          f"(F50 would need a transit shorter than {t_exp*1e3*f3['F50']/(E1.mean()*Pmax/(A_fly/100)):.0f} ms)")
