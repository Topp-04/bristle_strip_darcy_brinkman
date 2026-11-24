#!/usr/bin/env python3
"""Drag coefficient of the solid plate vs Janour (1951) and Blasius. Called by validation/plate_ReL120/Allrun.

Cd from forceCoeffs with Aref = 2 L Lz, i.e. Janour's Cx (drag per unit area of both sides).
Janour (1951) fit: Cx = 2.90 Re_L^-0.601 (10 < Re_L < 1000, +-3 %). Exit code 2 = still drifting.
"""
import glob, os, sys
import numpy as np
info = dict(l.split("=", 1) for l in open("case.info").read().split() if "=" in l)
ReL = float(info["ReL"])
files = sorted(glob.glob("postProcessing/forceCoeffs/*/*.dat"), key=os.path.getmtime)
if not files: sys.exit("no forceCoeffs output found in postProcessing/forceCoeffs/")
t, cd = [], []
for f in files:
    head = None
    for l in open(f):
        if l.startswith("#"):
            if "Cd" in l: head = l[1:].split()
            continue
        v = l.split()
        if head and len(v) == len(head):
            t.append(float(v[0])); cd.append(float(v[head.index("Cd")]))
t, cd = np.array(t), np.array(cd)
o = np.argsort(t); t, cd = t[o], cd[o]
last = cd[-1]
ref = cd[np.searchsorted(t, t[-1] - 500)] if t[-1] > 500 else np.nan
drift = abs(last - ref) / abs(last)
jan = 2.90 * ReL ** -0.601; bla = 1.328 / np.sqrt(ReL)
print(f"Re_L = {ReL:g}   iteration {t[-1]:.0f}:  Cd = {last:.5f}   drift {100 * drift:.3f} % per 500 it")
print(f"  Janour (1951) 2.90 Re^-0.601 = {jan:.5f}   ->  sim/Janour = {last / jan:.3f}  ({100 * (last / jan - 1):+.1f} %)")
print(f"  Blasius 1.328/sqrt(Re)       = {bla:.5f}   ->  sim/Blasius = {last / bla:.3f}")
if drift > 0.005:
    print("  -> still drifting"); sys.exit(2)
open("validation_result.csv", "w").write(f"ReL,Cd_sim,Cd_Janour,Cd_Blasius,drift_pct,iteration\n{ReL},{last:.6f},{jan:.6f},{bla:.6f},{100 * drift:.4f},{t[-1]:.0f}\n")
