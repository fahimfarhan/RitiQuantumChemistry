from warnings import simplefilter

import psi4

def dft_water():
    # ==============================================================
    # 1️⃣ Define the molecule (initial geometry, charge, spin)
    # ==============================================================
    # h2o = psi4.geometry("""
    #     0 1               # total charge = 0, multiplicity = 1 (singlet)
    #     O                 # oxygen atom at origin
    #     H 1 0.958         # first hydrogen 0.958 Å from O
    #     H 1 0.958 2 104.5 # second hydrogen 0.958 Å away, 104.5° angle
    #     units angstrom
    #     symmetry c1       # disable symmetry for generality
    #     """)

    h2o = psi4.geometry("""
        0 1
        O   0.000000   0.000000   0.000000
        H   0.000000   0.757000   0.586000
        H   0.000000  -0.757000   0.586000
        units angstrom
        symmetry c1
    """)

    # ==============================================================
    # 2️⃣ Choose the method and basis set
    # ==============================================================
    method = 'B3LYP'
    basis = '6-31G(d)'

    # ==============================================================
    # 3️⃣ Global computation options
    # ==============================================================
    psi4.set_options({
        'basis': basis,
        'd_convergence': 1e-8,
        'e_convergence': 1e-8,
        'guess': 'sad',  # initial guess: superposition of atomic densities
        'geom_maxiter': 50,  # maximum geometry optimization steps
        'print': 2,  # moderate verbosity

        # *** Recommended CPU Optimization ***
        'scf_type': 'df',  # Density Fitting is the biggest speedup for DFT/HF
        'freeze_core': True,  # Saves time by excluding core electrons from correlation
        # (only matters for post-HF like MP2, but good practice)
        'num_threads': -1  # Tells Psi4 to use ALL available CPU cores
    })

    # ==============================================================
    # 4️⃣ Single-point DFT energy (for reference)
    # ==============================================================
    E_sp = psi4.energy(method, molecule=h2o)
    print(f"\n🔹 Single-point DFT Energy ({method}/{basis}) = {E_sp:.10f} Hartree")

    # ==============================================================
    # 5️⃣ Geometry optimization (relax structure)
    # ==============================================================
    # Psi4 uses analytic gradients and the Berny algorithm.
    # This will iteratively update nuclear coordinates until forces vanish.
    print("\n🔹 Optimizing geometry ...")
    E_opt, wfn = psi4.optimize(method, molecule=h2o, return_wfn=True)
    print(f"✅ Optimized Energy = {E_opt:.10f} Hartree")

    # ==============================================================
    # 6️⃣ Property analysis (charges and dipole)
    # ==============================================================
    psi4.properties(method, properties=['MULLIKEN_CHARGES', 'dipole'], molecule=h2o)

    # ==============================================================
    # 7️⃣ Visualization output (for external viewers)
    # ==============================================================
    # Save MOLDEN file for visualization in programs like:
    #   - JMol, Molden, Avogadro, VMD
    # Contains: optimized geometry, orbitals, and electron density info.
    psi4.molden(wfn, 'water_optimized.molden')
    print("\n💾 Molden file saved as 'water_optimized.molden'")

    # ==============================================================
    # 8️⃣ Summary
    # ==============================================================
    print("""
        ================== Summary ==================
        Method/Basis  : {} / {}
        E(single pt)  : {:.10f} Hartree
        E(optimized)  : {:.10f} Hartree
        Molden file   : water_optimized.molden
        =============================================
        """.format(method, basis, E_sp, E_opt))

    # Define water molecule




    # Optimize geometry
    print("optimization start")
    opt_output = psi4.optimize(method, molecule=h2o)
    print(opt_output)
    print("optimization end")

    # Compute dipole moment and ESP charges
    # Run an energy calculation and get the wavefunction
    energy, wfn = psi4.energy(method, molecule=h2o, return_wfn=True)

    # Now compute dipole and ESP charges using the wavefunction
    psi4.oeprop(wfn, 'DIPOLE', 'ESP_CHARGES')

    # Write to Molden file
    psi4.molden(wfn, 'water.molden')
    pass


if __name__ == '__main__':
    dft_water()
    pass