import os
import psi4

# ---------- USER PATHS ----------
BASE = r"C:\Users\RITI\PycharmProjects\PythonProject1"
CELL = os.path.join(BASE, "cellobiose.xyz")
WATR = os.path.join(BASE, "water.xyz")
MEOH = os.path.join(BASE, "methanol.xyz")

# ---------- OUTPUT & SYSTEM SETTINGS ----------
psi4.core.set_output_file(os.path.join(BASE, "cellobiose_cluster.out"), False)
psi4.set_memory("8 GB")
psi4.core.set_num_threads(4)

# ---------- FUNCTION: READ .XYZ FILE ----------
def read_xyz(file):
    """Read an XYZ file and return a list of (atom, x, y, z)."""
    atoms = []
    with open(file, "r") as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]
        if len(lines) < 3:
            raise ValueError(f"❌ File {file} seems empty or invalid XYZ format.")
        try:
            int(lines[0])  # atom count
            start = 2
        except ValueError:
            start = 1
        for l in lines[start:]:
            parts = l.split()
            if len(parts) >= 4:
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    if not atoms:
        raise ValueError(f"❌ No atoms found in file {file}")
    return atoms

# ---------- LOAD FRAGMENTS ----------
cell = read_xyz(CELL)
wat = read_xyz(WATR)
meoh = read_xyz(MEOH)

# ---------- TRANSLATE TO SEPARATE FRAGMENTS ----------
def translate(coords, dx=0, dy=0, dz=0):
    return [(a, x + dx, y + dy, z + dz) for a, x, y, z in coords]

wat = translate(wat, dx=4.0)
meoh = translate(meoh, dx=-4.0)

cluster = cell + wat + meoh

# ---------- BUILD PSI4 MOLECULE STRING (NO HEADER) ----------
lines = ["0 1"]  # Charge 0, Multiplicity 1

for a, x, y, z in cluster:
    lines.append(f"{a:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")

lines.append("units angstrom")
lines.append("no_com")
lines.append("no_reorient")

mol_str = "\n".join(lines).strip() + "\n\n"

# ---------- DEBUG PRINT ----------
print("\n=== FINAL CLUSTER INPUT SENT TO PSI4 ===\n")
print(mol_str)

# ---------- CREATE PSI4 MOLECULE ----------
mol = psi4.geometry(mol_str)
print(f"✅ Cluster built successfully with {mol.natom()} atoms.\n")

# ---------- STAGE 1: HF/6-31G* PREOPTIMIZATION ----------
print("🔧 Stage 1: HF/6-31G* pre-optimization ...\n")
psi4.set_options({
    "basis": "6-31G*",
    "reference": "rhf",
    "scf_type": "pk",
    "diis": True,
    "soscf": True,
    "e_convergence": 1e-6,
    "d_convergence": 1e-6,
    "maxiter": 200,
    "opt_coordinates": "cartesian",
    "geom_maxiter": 200,
    "print": 2,
})

e1, wfn1 = psi4.optimize("hf", molecule=mol, return_wfn=True)
print(f"✅ Stage 1 energy (HF/6-31G*): {e1:.10f} Eh\n")

# ---------- STAGE 2: DFT OPTIMIZATION ----------
print("🚀 Stage 2: B3LYP-D3BJ/6-31++G(2d,2p) optimization ...\n")
psi4.set_options({
    "basis": "6-31++G(2d,2p)",
    "reference": "rhf",
    "scf_type": "pk",
    "diis": True,
    "soscf": True,
    "e_convergence": 1e-8,
    "d_convergence": 1e-8,
    "maxiter": 300,
    "opt_coordinates": "cartesian",
    "geom_maxiter": 200,
    "print": 2,
})

e2, wfn2 = psi4.optimize("b3lyp-d3bj", molecule=wfn1.molecule(), return_wfn=True)
print(f"✅ Stage 2 energy (B3LYP-D3BJ/6-31++G(2d,2p)): {e2:.10f} Eh\n")

# ---------- FREQUENCY ANALYSIS ----------
print("🔍 Harmonic frequency analysis ...\n")
psi4.set_options({
    "scf_type": "pk",
    "e_convergence": 1e-8,
    "d_convergence": 1e-8,
    "print": 2,
})
psi4.frequency("b3lyp-d3bj", molecule=wfn2.molecule())

# ---------- SAVE FINAL STRUCTURE ----------
xyz_out = os.path.join(BASE, "cluster_optimized.xyz")
wfn2.molecule().save_xyz_file(xyz_out, 12)
print(f"💾 Final optimized geometry saved to {xyz_out}\n")

# ---------- SUMMARY ----------
print("🎉 DFT Cluster Optimization Complete!")
print(f"HF Energy  (6-31G*):            {e1:.10f} Eh")
print(f"B3LYP-D3BJ/6-31++G(2d,2p):     {e2:.10f} Eh")
print("✅ Frequency analysis done.\n")