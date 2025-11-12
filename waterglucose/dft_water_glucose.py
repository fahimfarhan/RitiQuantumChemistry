import psi4

def dft_water_glucose():
    # ==============================================================
    # 1️⃣ Define the geometry of the water_glucose_molecule
    # ==============================================================

    # We'll use a small glucose fragment for simplicity (C6H12O6 is big).
    # You can replace this with a full geometry from PubChem or XYZ file later.

    # For demonstration, we’ll make a minimalistic model:
    # glucose oxygen interacting with water via H-bond

    water_glucose_molecule = psi4.geometry("""
    0 1
    C    0.0000    0.0000    0.0000
    O    1.2000    0.0000    0.0000
    H    1.5000    0.9000    0.0000
    H    1.5000   -0.9000    0.0000
    --
    O    2.7000    0.0000    0.0000
    H    3.3000    0.8000    0.0000
    H    3.3000   -0.8000    0.0000
    units angstrom
    symmetry c1
    """)

    # Define fragments (for BSSE correction)
    fragments = water_glucose_molecule.get_fragments()
    print(f"fragments = {fragments}")

    # ==============================================================
    # 2️⃣ Choose DFT method and basis set
    # ==============================================================

    method = 'B3LYP-D3BJ'  # includes empirical dispersion (important for H-bond)
    basis = '6-31G(d)'

    psi4.set_options({
        'basis': basis,
        'scf_type': 'df',
        'd_convergence': 1e-8,
        'e_convergence': 1e-8,
        'guess': 'sad'
    })

    # ==============================================================
    # 3️⃣ Compute total energy of water_glucose_molecule (with counterpoise correction)
    # ==============================================================

    E_complex = psi4.energy(f"{method}/{basis}", molecule=water_glucose_molecule, bsse_type='cp')

    # ==============================================================
    # 4️⃣ Compute energies of isolated fragments
    # ==============================================================

    # These are automatically available from counterpoise decomposition:
    # E_water = psi4.variable(fragments[1]) # ("FRAG: 1 CP ENERGY")
    # E_glucose = psi4.variable(fragments[0]) # ("FRAG: 0 CP ENERGY")
    E_glucose = psi4.variable("FRAG_0_CP")
    E_water = psi4.variable("FRAG_1_CP")
    # ==============================================================
    # 5️⃣ Compute interaction energy
    # ==============================================================

    E_interaction = E_complex - (E_glucose + E_water)

    # ==============================================================
    # 6️⃣ Print results
    # ==============================================================

    print("\n===== RESULTS =====")
    print(f"Total water_glucose_molecule energy (corrected): {E_complex:.10f} Hartree")
    print(f"Glucose fragment energy:          {E_glucose:.10f} Hartree")
    print(f"Water fragment energy:            {E_water:.10f} Hartree")
    print(f"Interaction energy (ΔE_int):      {E_interaction * 627.509:.3f} kcal/mol")

    # ==============================================================
    # 🧠 Mapping to physical steps
    # ==============================================================

    # 1. Define atomic coordinates → defines nuclear potential (V_ne)
    # 2. Compute electronic density ρ(r) using DFT (B3LYP-D3BJ)
    # 3. Include dispersion correction (D3BJ) for H-bond / van der Waals
    # 4. Counterpoise correction removes basis set overlap error
    # 5. Subtract monomer energies → gives binding strength

    psi4.optimize(method, molecule=water_glucose_molecule)

    penergy, wfn = psi4.energy(method, molecule=water_glucose_molecule, return_wfn=True)

    # Now compute dipole and ESP charges using the wavefunction
    psi4.oeprop(wfn, 'DIPOLE', 'ESP_CHARGES')

    # Write to Molden file
    psi4.molden(wfn, 'water_glucose.molden')
    pass

if __name__ == '__main__':
    dft_water_glucose()
    pass