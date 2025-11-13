import sys

from commoncode import *

h2o = """
O  0.94921  -0.06421  -0.05148
H  0.67036  -0.78825  -0.63196
H  1.91710  -0.10043  -0.08052
"""

methanol = """
O  0.70790  0.00000  0.00000
C -0.70790  0.00000  0.00000
H -1.07320 -0.76900  0.68520
H -1.07310 -0.19470 -1.01130
H -1.06320  0.97860  0.33120
H  0.99360 -0.88040 -0.29800
"""

cellobiose = """
C   -3.9780   -2.4560   0.0000
C   -2.6820   -1.7890   0.0000
C   -1.5210   -2.7460   0.0000
O   -1.4910   -4.1450   0.0000
C   -0.2210   -2.0070   0.0000
C   -0.2600   -0.5090   0.0000
O   -1.3660    0.2500   0.0000
C    1.0650    0.1880   0.0000
C    2.2250   -0.7690   0.0000
O    3.4900   -0.0470   0.0000
C    2.2680   -2.2750   0.0000
O    3.3980   -3.0450   0.0000
O    0.9340   -2.9400   0.0000
H   -4.8970   -1.8440   0.0000
H   -4.0390   -3.5340   0.0000
H   -2.7230   -0.7080   0.0000
H    0.7330   -0.0930   0.0000
H    1.1510    1.2700   0.0000
H    2.2060   -2.7930   0.9320
H    2.2060   -2.7930  -0.9320
H    3.8530   -3.5520   0.0000
H    0.8950   -3.5200   0.9320
H    0.8950   -3.5200  -0.9320
H   -0.1690    0.7070   0.0000
O   -5.1900   -0.9400   0.0000
O   -3.8720    0.1800   0.0000
C   -2.4860    0.9500   0.0000
C   -1.2750   -0.0550   0.0000
C    0.0400    0.7030   0.0000
O    1.1570   -0.0190   0.0000
C    2.3660    0.7140   0.0000
O    2.8860   -0.4960   0.0000
C    4.3060   -0.4960   0.0000
C    4.8260    0.7140   0.0000
O    4.3060    1.9240   0.0000
C    2.8860    1.9240   0.0000
O    2.3660    3.1340   0.0000
C    0.9460    3.1340   0.0000
O    0.4260    1.9240   0.0000
H   -6.0150   -1.4400   0.0000
H   -4.0180    1.2600   0.0000
H   -2.4780    2.0340   0.0000
H    5.8730    0.7140   0.0000
H    4.8260   -1.4150   0.0000
H    0.9460    4.0500   0.0000
H    2.3660    4.0500   0.0000
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
        {cellobiose}
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
    cellobiose_geom = psi4.geometry(f"""
        0 1
        {cellobiose}
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
    E_cellobiose_opt, _ = optimize_and_get_energy(cellobiose_geom, METHOD, BASIS, "cellobiose_Monomer")

    # 4. Property Analysis and Output for Trimer
    calculate_properties_and_save(wfn_trimer, "Cluster")

    # 5. Calculate Binding Energy (Total Interaction Energy)
    # E_binding = E(Trimer) - E_h2o_opt - E_methanol_opt - E_cellobiose_opt
    E_binding = E_trimer_opt - (E_h2o_opt + E_methanol_opt + E_cellobiose_opt)
    E_binding_kcal = E_binding * psi4.constants.hartree2kcalmol

    # 6. Final Summary
    summary_message = f"""
            ================== Cluster Final Summary ==================
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
