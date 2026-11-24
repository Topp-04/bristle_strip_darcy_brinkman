#!/usr/bin/env python3
"""Measure the drag of the porous strip (run from inside a strip case folder).

    python3 ../../../scripts/strip_cf.py              measure the latest solution, append to calib_log.csv
    python3 ../../../scripts/strip_cf.py --set-da X   write Da = X into fvModels (before the next run)

Interior C_F = drag per unit length of the strip over cylinders 25-74, divided by 1/2 rho U^2:
    C_F = 2 (nu/K) * (integral of u over those strip cells) / (window length * Lz * U^2)
Exit code 2 = solution still changing (C_F drift > 0.5 % per write), so Allrun runs more iterations.
"""
import os, re, sys, subprocess
import numpy as np

D = 0.01; LZ = D
info = dict(l.split("=", 1) for l in open("case.info").read().split() if "=" in l)
GD = int(info["GD"]); U = float(info["U"]); NU = float(info["nu"]); TARGET = float(info["CF_target"])
L = 100 * (GD + 1) * D + D                     # strip length

def set_da(da):                                # Darcy term coefficient in fvModels: implicit = -nu/K
    for f in ["constant/fvModels", "system/fvModels"]:
        txt = re.sub(r"(implicit\s+)[-+0-9.eE]+", lambda m: f"{m.group(1)}{-NU / (da * D**2):.8e}", open(f).read())
        open(f, "w").write(txt)

def K_used(t):                                 # K of the run that produced time t (snapshot saved by Allrun)
    snaps = sorted(int(f.split(".")[-1]) for f in os.listdir(".") if f.startswith("fvModels.used."))
    snap = min(s for s in snaps if s >= int(t))
    coeff = float(re.search(r"implicit\s+([-+0-9.eE]+)", open(f"fvModels.used.{snap}").read()).group(1))
    return -NU / coeff

def read_field(path, ncomp):                   # internalField of an ASCII OpenFOAM field
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", txt)
    a = np.array(txt[m.end():txt.index("\n)", m.end())].replace("(", " ").replace(")", " ").split(), float)
    return a.reshape(-1, ncomp) if ncomp > 1 else a

def cf_at(t, K):
    if not os.path.exists("cells.npz"):        # cell centres C and volumes V (written once)
        for fn in ["writeCellCentres", "writeCellVolumes"]:
            subprocess.run(f"foamPostProcess -func {fn} -time {t} > log.{fn} 2>&1", shell=True, check=True)
        np.savez("cells.npz", C=read_field(f"{t}/C", 3), V=read_field(f"{t}/V", 1))
    z = np.load("cells.npz"); C, V = z["C"], z["V"]
    u = read_field(f"{t}/U", 3)[:, 0]
    strip = (np.abs(C[:, 1]) < D / 2) & (C[:, 0] > 0) & (C[:, 0] < L)
    p = (GD + 1) * D
    win = strip & (C[:, 0] > D / 2 + 25 * p) & (C[:, 0] < D / 2 + 75 * p)        # cylinders 25-74
    cf = 2 * (NU / K) * np.sum(u[win] * V[win]) / (LZ * 50 * p) / U**2
    # C_F along the strip, column by column
    xs, inv = np.unique(np.round(C[strip, 0], 9), return_inverse=True)
    ubar = np.bincount(inv, weights=u[strip] * V[strip]) / np.bincount(inv, weights=V[strip])
    np.savetxt("CF_x_strip.csv", np.c_[xs / D, 2 * (NU / K) * ubar * D / U**2, ubar / U], delimiter=",",
               header="x_over_D_from_strip_start,CF_strip,ubar_star_strip", comments="")
    return cf, np.mean(u[win]) / U

if "--set-da" in sys.argv:
    set_da(float(sys.argv[sys.argv.index("--set-da") + 1])); sys.exit(0)

ts = sorted((d for d in os.listdir(".") if re.fullmatch(r"[0-9]+", d) and d != "0"), key=int)
K = K_used(ts[-1])
cf_prev = cf_at(ts[-2], K)[0] if len(ts) > 1 else np.nan
cf, ubar = cf_at(ts[-1], K)
drift = abs(cf - cf_prev) / cf
err = (cf - TARGET) / TARGET
print(f"iteration {ts[-1]}: Da = {K / D**2:.6g}  C_F strip = {cf:.5f}  C_F resolved = {TARGET:.5f}  "
      f"error = {100 * err:+.1f} %  drift = {100 * drift:.3f} %")
if drift > 0.005:
    sys.exit(2)
new = not os.path.exists("calib_log.csv")
with open("calib_log.csv", "a") as f:
    if new: f.write("time,Da,K_m2,CF_strip,CF_target,error_pct,drift_pct,ubar_star_strip\n")
    f.write(f"{ts[-1]},{K / D**2:.6g},{K:.6g},{cf:.6f},{TARGET:.6f},{100 * err:.3f},{100 * drift:.3f},{ubar:.5f}\n")
