import psi4


def dft_water_dimer():
    psi4.set_memory('1 GB')  # Increased memory for the larger system

    # 1. Define the Water Dimer Geometry (using Z-matrix for clarity)
    # This geometry is a classic H-bonded dimer configuration
    dimer_geom = psi4.geometry("""
        0 1                 # Charge 0, Multiplicity 1 (Molecule A)
        O
        H 1 0.957
        H 1 0.957 2 104.5
        --
        0 1                 # Charge 0, Multiplicity 1 (Molecule B)
        O 1 2.98 2 180.0 3 0.0  # O of B is 2.98 Å from O of A
        H 4 0.957 5 104.5 2 180.0  # H on O-H donor axis
        H 4 0.957 5 104.5 2 60.0   # The other H
    """)

    # 2. Choose the method and basis set
    method = 'B3LYP-D3'  # **Crucial Upgrade**: Use a dispersion-corrected DFT
    basis = '6-31G(d)'

    psi4.set_options({
        'basis': basis,
        'd_convergence': 1e-8,
        'e_convergence': 1e-8,
        'guess': 'sad',
        'geom_maxiter': 100,  # More steps needed for a complex

        # *** Recommended CPU Optimization ***
        'scf_type': 'df',  # Density Fitting is the biggest speedup for DFT/HF
        'freeze_core': True,  # Saves time by excluding core electrons from correlation
        # (only matters for post-HF like MP2, but good practice)
        'num_threads': -1  # Tells Psi4 to use ALL available CPU cores
    })

    # --- 3. Geometry Optimization of the Complex (A...B) ---
    print("\n🔬 Optimizing Water Dimer (A...B) geometry...")
    E_dimer_opt, wfn_dimer = psi4.optimize(f'{method}/{basis}', molecule=dimer_geom, return_wfn=True)
    print(f"✅ Optimized Dimer Energy = {E_dimer_opt:.8f} Hartree")

    # --- 4. Single-Point Energy of the Monomers (A and B) ---
    # The complex optimization also calculates the energy of the monomers
    # at their isolated, optimized geometry. This is the simplest way.

    # Alternatively, you could extract the geometry of A and B from the
    # optimized dimer (wfn_dimer.molecule().pyrian_geometry()) and optimize them separately.
    # For simplicity, we'll run a single-point on one monomer.

    # We must calculate the energy of the *isolated monomer*
    monomer_geom = psi4.geometry("""
        0 1
        O
        H 1 0.957
        H 1 0.957 2 104.5
    """)

    # Optimize the single water monomer
    print("\n🔬 Optimizing Water Monomer (A) geometry...")
    E_monomer_opt = psi4.optimize(f'{method}/{basis}', molecule=monomer_geom)
    print(f"✅ Optimized Monomer Energy = {E_monomer_opt:.8f} Hartree")

    # --- 5. Calculate Interaction Energy (Binding Energy) ---
    # E(A) and E(B) are identical due to symmetry, so E_A + E_B = 2 * E_monomer_opt
    E_binding = E_dimer_opt - (2 * E_monomer_opt)

    # Convert to a common chemistry unit (kcal/mol) for easier comparison
    E_binding_kcal = E_binding * psi4.constants.hartree2kcalmol

    print(f"""
        ================== Dimer Summary ==================
        Method/Basis  : {method} / {basis}
        E(Dimer)      : {E_dimer_opt:.8f} Hartree
        E(Monomer)    : {E_monomer_opt:.8f} Hartree

        🔹 Interaction Energy (2E_monomer - E_dimer) 🔹
        Interaction (Hartree) : {E_binding:.6f} Hartree
        Interaction (kcal/mol): {E_binding_kcal:.2f} kcal/mol
        =================================================
    """)

def dft_water_dimer_optimized():
    # --- 0. Setup and CPU Optimization ---

    # Set memory for the calculation (1 GB should be plenty for a dimer)
    psi4.set_memory('1 GB')

    # Set number of threads to use all available CPU cores
    # This is key for your multi-core CPU (Asus G14)
    psi4.set_num_threads(4)

    # Define method and basis set
    method = 'B3LYP-D3'  # DFT with Grimme's dispersion correction
    basis = '6-31G(d)'

    psi4.set_options({
        'basis': basis,
        'd_convergence': 1e-8,
        'e_convergence': 1e-8,
        'guess': 'sad',
        'geom_maxiter': 50, # use 100 or more in final
        'print': 2,

        # *** CPU Optimization ***
        'scf_type': 'df',  # Density Fitting is the biggest speedup for DFT/HF
        'freeze_core': True,  # Saves time by excluding core electrons from correlation (only matters for post-HF like MP2, but good practice)
    })

    # --- 1. Define the Water Dimer Geometry ---
    # Define the complex with a geometry that favors H-bonding (Molecule A -- Molecule B)
    # The '--' separates fragments for methods like SAPT, but is good practice here.

    """
            0 1                 # Charge 0, Multiplicity 1 (Fragment 1)
        O
        H 1 0.957
        H 1 0.957 2 104.5
        --
        0 1                 # Charge 0, Multiplicity 1 (Fragment 2)
        O 1 2.98 2 180.0 3 0.0  # O of B is 2.98 Å from O of A
        H 4 0.957 5 104.5 2 180.0
        H 4 0.957 5 104.5 2 60.0
    """

    dimer_geom = psi4.geometry("""
        0 1
        O                 
        H 1 0.957         
        H 1 0.957 2 104.5
        --
        0 1
        O 1 2.98 2 180.0 3 0.0
        H 4 0.957 1 104.5 2 180.0
        H 4 0.957 5 104.5 1 60.0
    """)

    # --- 2. Geometry Optimization of the Complex (A...B) ---
    print("\n🔬 Optimizing Water Dimer (A...B) geometry...")
    E_dimer_opt, wfn_dimer = psi4.optimize(f'{method}/{basis}', molecule=dimer_geom, return_wfn=True)
    print(f"✅ Optimized Dimer Energy = {E_dimer_opt:.8f} Hartree")

    # --- 3. Property Analysis and Visualization Output ---

    # Use the optimized wavefunction (wfn_dimer) to calculate properties
    print("\nCalculating properties (Dipole and Charges)...")

    # Calculate properties: Dipole Moment and Electrostatic Potential (ESP) charges
    psi4.oeprop(wfn_dimer, 'DIPOLE', 'ESP_CHARGES')

    # Generate MOLDEN file for visualization
    molden_file = 'water_dimer_optimized.molden'
    psi4.molden(wfn_dimer, molden_file)
    print(f"\n💾 Molden file saved as '{molden_file}'")

    # --- 4. Monomer Optimization and Binding Energy ---

    # Optimize a single water monomer (Fragment A)
    monomer_geom = psi4.geometry("""
        0 1
        O
        H 1 0.957
        H 1 0.957 2 104.5
    """)
    print("\n🔬 Optimizing Water Monomer (A) geometry...")
    E_monomer_opt = psi4.optimize(f'{method}/{basis}', molecule=monomer_geom)

    # Calculate binding energy (subtract 2x monomer energy from dimer energy)
    E_binding = E_dimer_opt - (2 * E_monomer_opt)
    E_binding_kcal = E_binding * psi4.constants.hartree2kcalmol

    # --- 5. Final Summary ---
    print(f"""
        ================== Water Dimer Summary ==================
        Method/Basis  : {method} / {basis} (CPU Optimized)
        E(Optimized Dimer) : {E_dimer_opt:.8f} Hartree
        E(Optimized Monomer): {E_monomer_opt:.8f} Hartree

        🔹 Interaction Energy (Binding Energy) 🔹
        Interaction (Hartree) : {E_binding:.6f} Hartree
        Interaction (kcal/mol): {E_binding_kcal:.2f} kcal/mol
        =======================================================
    """)


if __name__ == '__main__':
    # Set the output file name for detailed log (optional, but helpful)
    # dft_water_dimer()

    psi4.set_output_file('water_dimer.log', False)
    dft_water_dimer_optimized()

