import sys

from commoncode import *

h2o = """
O            4.863548169202     0.135042546524     0.000003358471
H            4.115055602571    -0.475514748018     0.000001289630
H            5.642405360243    -0.435209389734     0.000004323706
"""

methanol = """
C            0.716615969945     0.072467437735    -0.000004833050
H            0.579388382383     0.773185355040     1.036950473417
H            0.201417676647    -0.277022476256    -0.000006255870
H            0.579392564603     0.773187391365    -1.036959312668
O            1.815277371461    -0.551748051113    -0.000003293956
H            1.886865400471    -1.511408686766    -0.000004098699
"""

co2 = """
C           -2.922436901136     0.159084015735    -0.000011621813
O           -4.087858826501     0.099718874015     0.000022974848
O           -1.755478330631     0.215903607308    -0.000009838598
"""

def define_complex_molecule_geometry():
    """
    Defines the cyclic water trimer (H2O)3 complex (9 atoms)
    using correct atom indexing for the Z-matrix format.
    """

    three_molecules = f"""
        0 1
        {h2o}
        --
        {methanol}
        --
        {co2}
        """

    print(three_molecules)

    trimer_geom = psi4.geometry(three_molecules)
    return trimer_geom

def three_molecule_analysis():
    # 1. Setup Environment
    setup_psi4_environment(METHOD, BASIS, THREADS)

    trimer_geom = define_complex_molecule_geometry()

    # Monomer geometry definition (used for isolated optimization). hmm. looks like monomer can have a different geometry, no need to share with the trimer above
    h2o_geom = psi4.geometry(f"""
        0 1
        {h2o}
    """)
    methanol_geom = psi4.geometry(f"""
        0 1
        {methanol}
    """)
    co2_geom = psi4.geometry(f"""
        0 1
        {co2}
    """)

    # 3. Perform Calculations
    E_trimer_opt = 0.0
    E_monomer_opt = 0.0

    # Trimer Optimization
    E_trimer_opt, wfn_trimer = optimize_and_get_energy(trimer_geom, METHOD, BASIS, "Three_Molecules")

    # CRITICAL CHECK: If trimer fails, stop the analysis.
    if E_trimer_opt == INFINITY:
        my_log_print("\n\n!!! CRITICAL FAILURE: Trimer optimization failed. Terminating analysis. !!!")
        sys.exit(1)

    # Monomer Optimization (Only need to do one, as they are identical)
    E_h2o_opt, _ = optimize_and_get_energy(h2o_geom, METHOD, BASIS, "Water_Monomer")
    E_methanol_opt, _ = optimize_and_get_energy(methanol_geom, METHOD, BASIS, "Methanol_Monomer")
    E_co2_opt, _ = optimize_and_get_energy(co2_geom, METHOD, BASIS, "CO2_Monomer")

    # 4. Property Analysis and Output for Trimer
    calculate_properties_and_save(wfn_trimer, "Three_Molecules")

    # 5. Calculate Binding Energy (Total Interaction Energy)
    # E_binding = E(Trimer) - E_h2o_opt - E_methanol_opt - E_co2_opt
    E_binding = E_trimer_opt - (E_h2o_opt + E_methanol_opt + E_co2_opt)
    E_binding_kcal = E_binding * psi4.constants.hartree2kcalmol

    # 6. Final Summary
    summary_message = f"""
            ================== Three Molecule Final Summary ==================
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
        psi4.set_output_file('three_molecule_analysis.log', False)
        three_molecule_analysis()
    except Exception as x:
        my_log_print(str(x))
    finally:
        psi4.core.clean()
        end = datetime.now()
        delta = end - start
        my_log_print(f"total run time {delta}.")
        my_log_print(f"threads count = {THREADS}")
    pass
