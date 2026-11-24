#!/usr/bin/env python3
"""Generate the six porous-strip cases (G/D = 2, 10; Re = 0.4, 0.5, 0.6) in cases/strip/.

The cases are already in the repository; run this only to regenerate them or to make new ones.
    python3 scripts/make_strip_cases.py            (from the repository root)
Mesh: blockMesh, 15 blocks; the middle block is the porous strip (cellZone 'strip').
Domain as in the resolved thesis simulations: 200D from the end cylinders, 200D to top/bottom.
"""
import os, math

D = 0.01                      # bristle diameter [m]
LZ = D                        # one-cell-thick 2D domain
NU = 6.3425e-4                # glycerin, kinematic viscosity [m^2/s]  (Re 0.4 <-> U = 0.02537 m/s)
NP = 8                        # MPI ranks (Allrun can override)
DX_END = 0.1 * D              # cell length at the strip ends
DX_MID = 1.0 * D              # cell length in the middle of the strip
DX_FAR = 4.0 * D              # cell length at inlet/outlet
DY_NEAR_OUT = 0.5 * D         # cell height 10D away from the strip
DY_FAR = 8.0 * D              # cell height at top/bottom boundaries
H_NEAR = 10.0 * D             # extent of the near-strip refinement band
WRITE_EVERY = 250             # iterations between writes (convergence is checked over this)

def read_geom(gd):
    """Upstream distance and half-height of the domain (as in the thesis): 200D each."""
    return 200 * D, 200 * D

# (G/D, Re): (plateau Cd of the resolved row, first-guess Da0)  -- see data/reference/Da0_strip_cases.csv
DATA = {(2, 0.4): (2.0779, 0.05193), (2, 0.5): (1.6765, 0.05204), (2, 0.6): (1.4062, 0.05216),
        (10, 0.4): (8.8584, 1.006), (10, 0.5): (7.0973, 1.005), (10, 0.6): (5.9220, 1.003)}

def HDR(obj, loc, cls="dictionary"):
    return ("FoamFile\n{\n    format      ascii;\n    class       %s;\n    location    \"%s\";\n"
            "    object      %s;\n}\n\n" % (cls, loc, obj))

def write(path, txt):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write(txt)

def ncells(L, a, b):
    """Number of cells for a geometric series of total length L from size a to size b."""
    if abs(a - b) < 1e-12 * L:
        return max(1, math.ceil(L / a))
    q = (L - a) / (L - b)
    return max(1, math.ceil(1 + math.log(b / a) / math.log(q)))

def block_mesh_dict(gd, da_mesh):
    up, H = read_geom(gd)
    L = 100 * (gd + 1) * D + D                     # strip length: first to last cylinder edge
    x_in, x_out = -(up - D / 2), L + (up - D / 2)  # inlet/outlet 'up' from the end cylinder centres
    ny_s = int(min(250, max(20, math.ceil(10 / math.sqrt(da_mesh)))))   # >=10 cells per sqrt(K)
    dy_s = D / ny_s
    X = [x_in, 0.0, L, x_out]
    Y = [-H, -H_NEAR, -D / 2, D / 2, H_NEAR, H]
    nx = [ncells(-x_in, DX_END, DX_FAR), 2 * ncells(L / 2, DX_END, DX_MID), ncells(x_out - L, DX_END, DX_FAR)]
    ny_near = ncells(H_NEAR - D / 2, dy_s, DY_NEAR_OUT)
    ny_far = ncells(H - H_NEAR, DY_NEAR_OUT, DY_FAR)
    ny = [ny_far, ny_near, ny_s, ny_near, ny_far]
    gx = [f"{DX_END / DX_FAR:.6g}",
          f"((0.5 0.5 {DX_MID / DX_END:.6g}) (0.5 0.5 {DX_END / DX_MID:.6g}))",
          f"{DX_FAR / DX_END:.6g}"]
    gy = [f"{DY_NEAR_OUT / DY_FAR:.6g}", f"{dy_s / DY_NEAR_OUT:.6g}", "1",
          f"{DY_NEAR_OUT / dy_s:.6g}", f"{DY_FAR / DY_NEAR_OUT:.6g}"]
    nxv, nyv = len(X), len(Y)
    vid = lambda i, j, k: k * nxv * nyv + j * nxv + i
    verts = [f"    ({X[i]:.8g} {Y[j]:.8g} {k * LZ:.8g})" for k in (0, 1) for j in range(nyv) for i in range(nxv)]
    blocks = []
    for j in range(nyv - 1):
        for i in range(nxv - 1):
            h = " ".join(str(v) for v in (vid(i, j, 0), vid(i + 1, j, 0), vid(i + 1, j + 1, 0), vid(i, j + 1, 0),
                                          vid(i, j, 1), vid(i + 1, j, 1), vid(i + 1, j + 1, 1), vid(i, j + 1, 1)))
            zone = " strip" if (i == 1 and j == 2) else ""
            blocks.append(f"    hex ({h}){zone} ({nx[i]} {ny[j]} 1) simpleGrading ({gx[i]} {gy[j]} 1)")
    face = lambda a, b, c, d: f"            ({a} {b} {c} {d})"
    inlet = [face(vid(0, j, 0), vid(0, j, 1), vid(0, j + 1, 1), vid(0, j + 1, 0)) for j in range(nyv - 1)]
    outlet = [face(vid(nxv - 1, j, 0), vid(nxv - 1, j + 1, 0), vid(nxv - 1, j + 1, 1), vid(nxv - 1, j, 1)) for j in range(nyv - 1)]
    sym = [face(vid(i, 0, 0), vid(i + 1, 0, 0), vid(i + 1, 0, 1), vid(i, 0, 1)) for i in range(nxv - 1)] + \
          [face(vid(i, nyv - 1, 0), vid(i, nyv - 1, 1), vid(i + 1, nyv - 1, 1), vid(i + 1, nyv - 1, 0)) for i in range(nxv - 1)]
    fb = [face(vid(i, j, 0), vid(i, j + 1, 0), vid(i + 1, j + 1, 0), vid(i + 1, j, 0)) for j in range(nyv - 1) for i in range(nxv - 1)] + \
         [face(vid(i, j, 1), vid(i + 1, j, 1), vid(i + 1, j + 1, 1), vid(i, j + 1, 1)) for j in range(nyv - 1) for i in range(nxv - 1)]
    patch = lambda name, typ, faces: f"    {name}\n    {{\n        type {typ};\n        faces\n        (\n" + "\n".join(faces) + "\n        );\n    }\n"
    total = sum(nx) * sum(ny)
    txt = (HDR("blockMeshDict", "system") +
           f"// porous strip G/D = {gd}: strip length {L / D:.0f}D, thickness D, {ny_s} cells across the strip\n"
           f"// cells: {sum(nx)} x {sum(ny)} = {total}\n\nconvertToMeters 1;\n\nvertices\n(\n" + "\n".join(verts) +
           "\n);\n\nblocks\n(\n" + "\n".join(blocks) + "\n);\n\nboundary\n(\n" +
           patch("inlet", "patch", inlet) + patch("outlet", "patch", outlet) +
           patch("sym", "symmetry", sym) + patch("frontAndBack", "empty", fb) + ");\n")
    info = dict(L=L, x_in=x_in, x_out=x_out, H=H, ny_strip=ny_s, cells=total)
    return txt, info

def fv_models(nu, K):
    # Darcy drag  -nu/K * U  in the strip only, added implicitly (kinematic form: eqn -= Sp(nu/K, U))
    return HDR("fvModels", "constant") + f"""// Darcy-Brinkman strip:  K = {K:.6e} m^2,  Da = K/D^2 = {K / D**2:.6e}
stripDarcy
{{
    type            semiImplicitSource;
    selectionMode   cellZone;      // OpenFOAM <= 11 keyword
    select          cellZone;      // newer keyword (unused entries are ignored)
    cellZone        strip;
    volumeMode      specific;
    sources
    {{
        U
        {{
            explicit    (0 0 0);
            implicit    {-nu / K:.8e};   // -nu/K [1/s]
        }}
    }}
}}
"""

def control_dict(end):
    return HDR("controlDict", "system") + f"""solver          incompressibleFluid;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {end};
deltaT          1;
writeControl    timeStep;
writeInterval   {WRITE_EVERY};
purgeWrite      3;
writeFormat     ascii;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable true;
"""

FV_SCHEMES = HDR("fvSchemes", "system") + """ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         Gauss linear;
    div(phi,U)      bounded Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
"""

FV_SOLUTION = HDR("fvSolution", "system") + """solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-10;
        relTol          0.05;
    }
    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-12;
        relTol          0.1;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;          // SIMPLEC: faster convergence for this steady viscous flow
}
relaxationFactors
{
    equations { U 0.9; ".*" 0.9; }
}
"""

def decompose_dict(n):
    return HDR("decomposeParDict", "system") + f"numberOfSubdomains {n};\nmethod          scotch;\n"

def field_U(U):
    return HDR("U", "0", "volVectorField") + f"""dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({U:.8g} 0 0);
boundaryField
{{
    inlet           {{ type fixedValue; value uniform ({U:.8g} 0 0); }}
    outlet          {{ type zeroGradient; }}
    sym             {{ type symmetry; }}
    frontAndBack    {{ type empty; }}
}}
"""

FIELD_P = HDR("p", "0", "volScalarField") + """dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    inlet           { type zeroGradient; }
    outlet          { type fixedValue; value uniform 0; }
    sym             { type symmetry; }
    frontAndBack    { type empty; }
}
"""

PHYS = HDR("physicalProperties", "constant") + f"viscosityModel  constant;\n\nnu              [0 2 -1 0 0 0 0] {NU:.6e};\n"
MT = HDR("momentumTransport", "constant") + "simulationType  laminar;\n"

if __name__ == "__main__":
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cases", "strip")
    for (gd, re_), (cd, da0) in DATA.items():
        case = os.path.join(root, f"GD{gd}_Re{re_:g}")
        U = re_ * NU / D
        bmd, info = block_mesh_dict(gd, da0 / 4)        # fine enough for K down to K0/4
        files = {"system/blockMeshDict": bmd, "system/controlDict": control_dict(4000),
                 "system/fvSchemes": FV_SCHEMES, "system/fvSolution": FV_SOLUTION,
                 "system/decomposeParDict": decompose_dict(NP),
                 "constant/physicalProperties": PHYS, "constant/momentumTransport": MT,
                 "constant/fvModels": fv_models(NU, da0 * D**2), "system/fvModels": fv_models(NU, da0 * D**2),
                 "0/U": field_U(U), "0/p": FIELD_P,
                 "case.info": f"GD={gd}\nRe={re_}\nU={U:.8g}\nnu={NU:.6e}\nCd_plateau={cd}\nCF_target={cd / (1 + gd):.6f}\nDa0={da0}\n"}
        for p, t in files.items():
            write(os.path.join(case, p), t)
        print(f"{case}: U = {U:.5f} m/s, Da0 = {da0}, target C_F = {cd / (1 + gd):.4f}, {info['cells']} cells")
