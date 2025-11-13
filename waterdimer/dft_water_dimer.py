from commoncode import *
# --- 3. CORE CALCULATION FUNCTIONS ---

def define_dimer_geometry():
    """Defines the water dimer Z-matrix and returns the psi4.Molecule object."""
    # Corrected Z-matrix for a hydrogen-bonded water dimer complex
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
    return dimer_geom




# --- 4. MAIN EXECUTION FUNCTION ---

def run_water_dimer_analysis():
    """Executes the full quantum chemistry workflow."""

    # 1. Setup Environment
    setup_psi4_environment(METHOD, BASIS, THREADS)

    # 2. Define Geometries
    dimer_geom = define_dimer_geometry()

    # Monomer geometry definition (used for isolated optimization)
    monomer_geom = psi4.geometry("""
        0 1
        O
        H 1 0.957
        H 1 0.957 2 104.5
    """)

    # 3. Perform Calculations
    E_dimer_opt = 0.0
    wfn_dimer = None

    E_monomer_opt = 0.0
    wfn_monomer = None

    try:
        E_dimer_opt, wfn_dimer = optimize_and_get_energy(dimer_geom, METHOD, BASIS, "Dimer")
    except Exception as x:
        my_log_print(x)
    try:
        E_monomer_opt, wfn_monomer = optimize_and_get_energy(monomer_geom, METHOD, BASIS, "Monomer")
    except Exception as x:
        my_log_print(x)

    # 4. Property Analysis and Output for Dimer
    try:
        calculate_properties_and_save(wfn_dimer, "Dimer")
    except Exception as x:
        my_log_print(x)
    # 5. Calculate Binding Energy
    E_binding = E_dimer_opt - (2 * E_monomer_opt)
    E_binding_kcal = E_binding * psi4.constants.hartree2kcalmol

    # 6. Final Summary
    my_log_print(f"""
        ================== Water Dimer Final Summary ==================
        Method/Basis  : {METHOD} / {BASIS} (CPU Optimized)
        E(Optimized Dimer) : {E_dimer_opt:.8f} Hartree
        E(Optimized Monomer): {E_monomer_opt:.8f} Hartree

        🔹 Interaction Energy (Binding Energy) 🔹
        Interaction (Hartree) : {E_binding:.6f} Hartree
        Interaction (kcal/mol): {E_binding_kcal:.2f} kcal/mol
        ===============================================================
    """)
    pass


if __name__ == '__main__':
    # Set the main log file name
    start = datetime.now()
    psi4.set_output_file('water_dimer_analysis.log', False)
    run_water_dimer_analysis()
    end = datetime.now()
    delta = end - start
    my_log_print(f"total run time {delta} seconds.")
    my_log_print(f"threads count = {THREADS}")
    pass
