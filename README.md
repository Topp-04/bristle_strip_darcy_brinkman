# Porous-strip (Darcy–Brinkman) model of flow along a bristled-wing array

Term project, *Modeling with Open Source Software*, Skoltech, Fall 2026 — Top D. T. Koryom.

A row of 101 bristles (cylinders, diameter `D = 0.01 m`, gap `G`) with the flow running along the row is
replaced by a porous strip of thickness `D`. Inside the strip a Darcy drag `−(ν/K)u` is added to the
Navier–Stokes equations. The strip's drag per unit length, `C_F = C_d D/(G+D)`, is compared with resolved
simulations of the real row (from my MSc thesis) for several values of the permeability `K`.
The setup is validated on a solid flat plate against the measurements of Janour (1951).

**Result:** whatever `K` is used, the strip's drag stays at the value of a solid plate, 2.1–4.1 times
lower than the real row, so no permeability reproduces the drag.

## What is in the repository

```
cases/strip/GD2_Re0.4 ... GD10_Re0.6   six porous-strip cases (G/D = 2, 10; Re = 0.4, 0.5, 0.6)
validation/plate_ReL120               solid plate, Re_L = 120, compared with Janour (1951)
scripts/                              Python helpers (below)
slurm/run_case.job                    run one case on a Slurm cluster
data/reference/                       data from the resolved thesis simulations
results/                              the small result files and figures used in the report
Allrun, Allclean                      run / clean all cases
```

Every case folder is a normal OpenFOAM case (`0/`, `constant/`, `system/`) plus:
- `Allrun`, `Allclean` — run and clean the case
- `case.info` — G/D, Re, inlet velocity U, viscosity, and the target `C_F` of the real row
- `Da_sequence` (strip cases only) — the Darcy numbers `Da = K/D²` that were run, in order

Scripts:
- `strip_cf.py` — measures the strip's `C_F` after a run and writes `calib_log.csv`; also sets a new K
- `plate_cd.py` — plate drag compared with Janour and Blasius
- `make_strip_cases.py` — regenerates the strip cases (all mesh settings are in this file)
- `plot_report_figures.py` — makes the report figures from `results/`

## Software

- **OpenFOAM 14** (openfoam.org), solver `foamRun` with `incompressibleFluid`; nothing is compiled
  (the Darcy term uses the standard `fvModel` `semiImplicitSource`)
- MPI (comes with OpenFOAM), Python 3 with `numpy` and `matplotlib`, bash
- Run on the Skoltech cluster with `module load openfoam/14` and 8 cores per case

## How to run

Load OpenFOAM first (`module load openfoam/14`, or `source /opt/openfoam14/etc/bashrc`).

**One strip case**
```bash
cd cases/strip/GD2_Re0.4
./Allrun                    # NP=4 ./Allrun to use 4 cores, NP=1 for serial
```
What `Allrun` does:
1. builds the mesh with `blockMesh` and splits it over the cores with `decomposePar`;
2. runs `foamRun` for the first Darcy number in `Da_sequence` (4000 iterations);
3. runs `strip_cf.py`, which measures `C_F` over cylinders 25–74 and adds a line to `calib_log.csv`;
   if `C_F` is still changing (more than 0.5 % between saves), 1500 more iterations are run first;
4. sets the next Darcy number and continues from the previous solution (1500 iterations), and so on.

**Validation case**
```bash
cd validation/plate_ReL120 && ./Allrun
```
Prints the plate drag next to Janour's fit `2.90 Re_L^-0.601` and Blasius `1.328 Re_L^-1/2`.

**All cases** — `./Allrun` in the top folder (one after another), or on the cluster, for each case:
`sbatch -J GD2_Re0.4 slurm/run_case.job cases/strip/GD2_Re0.4`

**Figures** — copy each case's `calib_log.csv` and `CF_x_strip.csv` into
`results/strip_calibration/<case>/`, then run `python3 scripts/plot_report_figures.py`.

**Clean** — `./Allclean` (in the top folder for all cases, or inside one case).

### Steps that are not automated
- The Darcy numbers in `Da_sequence` were chosen during the project: the first one (`Da0`) is an
  estimate from the thesis results (`data/reference/Da0_strip_cases.csv`); each next one came from
  comparing the measured `C_F` with the target. Since `C_F` hardly changed, K was then lowered by about
  10× at a time. `Allrun` simply repeats these values.
- `Da0` was computed from velocity profiles of the thesis simulations, which are not included here.

## Expected results

| Case | strip `C_F` (all Da run) | real row `C_F` | row / strip |
|---|---|---|---|
| G/D = 2, Re = 0.4 / 0.5 / 0.6 | 0.294–0.298 / 0.252–0.256 / 0.223–0.227 | 0.693 / 0.559 / 0.469 | 2.35 / 2.21 / 2.10 |
| G/D = 10, Re = 0.4 / 0.5 / 0.6 | 0.192–0.199 / 0.162–0.167 / 0.142–0.146 | 0.805 / 0.645 / 0.538 | 4.05 / 3.88 / 3.71 |
| Plate, Re_L = 120 | Cd = 0.171 (Janour 0.163, +5 %) | | |

## Notes
- The domain, boundary conditions and fluid (glycerin, ν = 6.3425×10⁻⁴ m²/s) are the same as in the
  thesis: inlet/outlet 200D from the end cylinders, symmetry planes 200D from the row.
- Mesh sizes: 0.43 M cells (G/D = 2), 0.82 M cells (G/D = 10), 0.38 M cells (plate, large domain).
- `fvModels` is kept in both `constant/` and `system/` because OpenFOAM versions read it from different
  places; both copies are always identical.

## References
- Z. Janour, *Resistance of a plate in parallel flow at low Reynolds numbers*, NACA TM 1316, 1951.
- G. Zhdanko, D. Kolomenskiy, S. E. Farisenkov, *Analytical criteria for bristle spacing on microinsect
  wings to maintain high aerodynamic force*, Phys. Fluids 37, 101707 (2025).
