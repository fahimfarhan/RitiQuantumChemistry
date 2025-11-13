import sys

from commoncode import *


def define_water_trimer_geometry():
    """
    Defines the cyclic water trimer (H2O)3 complex (9 atoms)
    using correct atom indexing for the Z-matrix format.
    """

    trimer_geom = psi4.geometry("""
    0 1
    # === Water 1 (Atoms 1-3) ===
    O
    H 1 0.957
    H 1 0.957 2 104.5
    --
    # === Water 2 (Atoms 4-6) ===
    # O2 is 2.85 A from O1 (Atom 1), 120 deg from H11 (Atom 2), 0 deg dihedral from H12 (Atom 3)
    O 1 2.85 2 120.0 3 0.0
    # H21 bonded to O2 (Atom 4), 104.5 deg from O1 (Atom 1), 180 deg dihedral from H11 (Atom 2)
    H 4 0.957 1 104.5 2 180.0
    # H22 bonded to O2 (Atom 4), 104.5 deg from H21 (Atom 5), 60 deg dihedral from O1 (Atom 1)
    H 4 0.957 5 104.5 1 60.0
    --
    # === Water 3 (Atoms 7-9) ===
    # O3 is 2.85 A from O2 (Atom 4), 120 deg from H21 (Atom 5), 0 deg dihedral from O1 (Atom 1)
    O 4 2.85 5 120.0 1 0.0
    # H31 bonded to O3 (Atom 7), 104.5 deg from O2 (Atom 4), 180 deg dihedral from H21 (Atom 5)
    H 7 0.957 4 104.5 5 180.0
    # H32 bonded to O3 (Atom 7), 104.5 deg from H31 (Atom 8), 60 deg dihedral from O2 (Atom 4)
    H 7 0.957 8 104.5 4 60.0
    """)
    return trimer_geom

def water_trimer_analysis():
    # 1. Setup Environment
    setup_psi4_environment(METHOD, BASIS, THREADS)

    trimer_geom = define_water_trimer_geometry()

    # Monomer geometry definition (used for isolated optimization). hmm. so monomer can have a different geometry, no need to share with the trimer above
    monomer_geom = psi4.geometry("""
            0 1
            O
            H 1 0.957
            H 1 0.957 2 104.5
        """)

    # 3. Perform Calculations
    E_trimer_opt = 0.0
    E_monomer_opt = 0.0

    # Trimer Optimization
    E_trimer_opt, wfn_trimer = optimize_and_get_energy(trimer_geom, METHOD, BASIS, "Water_Trimer")

    # CRITICAL CHECK: If trimer fails, stop the analysis.
    if E_trimer_opt == INFINITY:
        my_log_print("\n\n!!! CRITICAL FAILURE: Trimer optimization failed. Terminating analysis. !!!")
        sys.exit(1)

    # Monomer Optimization (Only need to do one, as they are identical)
    E_monomer_opt, _ = optimize_and_get_energy(monomer_geom, METHOD, BASIS, "Water_Monomer")

    # 4. Property Analysis and Output for Trimer
    calculate_properties_and_save(wfn_trimer, "Water_Trimer")

    # 5. Calculate Binding Energy (Total Interaction Energy)
    # E_binding = E(Trimer) - 3 * E(Monomer)
    E_binding = E_trimer_opt - (3 * E_monomer_opt)
    E_binding_kcal = E_binding * psi4.constants.hartree2kcalmol

    # 6. Final Summary
    summary_message = f"""
            ================== Water Trimer Final Summary ==================
            Method/Basis  : {METHOD} / {BASIS} (9 Atoms)
            E(Optimized Trimer)  : {E_trimer_opt:.8f} Hartree
            E(Optimized Monomer) : {E_monomer_opt:.8f} Hartree

            🔹 Total Interaction Energy (Binding Energy) 🔹
            Interaction (Hartree) : {E_binding:.6f} Hartree
            Interaction (kcal/mol): {E_binding_kcal:.2f} kcal/mol
            ================================================================
        """
    my_log_print(summary_message)

    pass

if __name__ == '__main__':
    # Set the main log file name
    start = datetime.now()
    try:
        psi4.set_output_file('water_trimer_analysis.log', False)
        water_trimer_analysis()
    except Exception as x:
        my_log_print(x)
    finally:
        psi4.core.clean()
        end = datetime.now()
        delta = end - start
        my_log_print(f"total run time {delta}.")
        my_log_print(f"threads count = {THREADS}")
    pass
