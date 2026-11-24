#!/usr/bin/env python3
"""Make the two result figures of the report from the stored results (no OpenFOAM needed).

    python3 scripts/plot_report_figures.py            (from the repository root)
Reads  results/strip_calibration/<case>/{calib_log,CF_x_strip}.csv   (copy fresh ones there after a rerun)
       data/reference/cd_all_cylinders.csv                            (resolved per-cylinder Cd, thesis)
Writes results/figures/fig_calibration.png, results/figures/fig_cf_along_row.png
"""
import csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RES = os.path.join(ROOT, "results", "strip_calibration")
OUT = os.path.join(ROOT, "results", "figures"); os.makedirs(OUT, exist_ok=True)

# Fig. 1: interior C_F of the strip vs Da, with the resolved-row target
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), dpi=200)
for ax, gd in zip(axs, [2, 10]):
    for i, re_ in enumerate([0.4, 0.5, 0.6]):
        rows = [r for r in csv.DictReader(open(os.path.join(RES, f"GD{gd}_Re{re_:g}", "calib_log.csv")))
                if float(r["drift_pct"]) < 0.5]                       # converged runs only
        ax.plot([float(r["Da"]) for r in rows], [float(r["CF_strip"]) for r in rows], "-o", color=f"C{i}",
                ms=5, lw=1.2, label=f"strip, Re = {re_}")
        ax.axhline(float(rows[0]["CF_target"]), color=f"C{i}", ls="--", lw=1, label=f"resolved row, Re = {re_}")
    ax.set_xscale("log"); ax.set_xlabel(r"Darcy number $\mathrm{Da}=K/D^2$"); ax.set_ylabel(r"interior $C_F$")
    ax.set_title(f"G/D = {gd}"); ax.grid(alpha=.3, which="both"); ax.set_ylim(0, None)
axs[1].legend(fontsize=7, loc="center left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_calibration.png"))

# Fig. 2: C_F along the row at Re = 0.4, resolved cylinders vs strip (smallest K), strip averaged per pitch
res = {}
for r in csv.DictReader(open(os.path.join(ROOT, "data", "reference", "cd_all_cylinders.csv"))):
    res.setdefault((int(float(r["GD"])), round(float(r["Re"]), 2)), []).append((int(r["cyl"]), float(r["Cd"])))

def pitch_avg(x, cf, p, n=101):
    o = np.argsort(x); x, cf = x[o], cf[o]; w = np.gradient(x)
    xs, vs = [], []
    for i in range(n):
        m = (x >= i * p) & (x < ((i + 1) * p if i < n - 1 else 1e9))
        if m.any(): xs.append(0.5 + i * p); vs.append(np.sum(cf[m] * w[m]) / np.sum(w[m]))
    return np.array(xs), np.array(vs)

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), dpi=200)
for ax, gd in zip(axs, [2, 10]):
    p = gd + 1
    a = np.loadtxt(os.path.join(RES, f"GD{gd}_Re0.4", "CF_x_strip.csv"), delimiter=",", skiprows=1)
    xs, vs = pitch_avg(a[:, 0], a[:, 1], p)
    cyl = sorted(res[(gd, 0.4)])
    ax.plot([0.5 + i * p for i, _ in cyl], [c / p for _, c in cyl], "o", ms=3, color="C0", label="resolved row (per cylinder)")
    ax.plot(xs, vs, "s", ms=3, color="C3", label="porous strip, solid limit (per pitch)")
    ax.axvspan(0.5 + 25 * p, 0.5 + 74 * p, color="0.9", zorder=0, label="interior window (cyl. 25–74)")
    ax.set_yscale("log"); ax.set_ylim(0.1, 5); ax.set_xlabel("x / D from the leading edge"); ax.set_ylabel(r"$C_F$")
    ax.set_title(f"G/D = {gd}, Re = 0.4"); ax.grid(alpha=.3, which="both")
axs[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_cf_along_row.png"))
print("wrote", OUT)
