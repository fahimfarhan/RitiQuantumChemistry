from commoncode import *

# --- 2. MOLECULE DEFINITIONS  ---


def define_molecules_for_test():
    """Defines a dictionary of small molecules for optimization and analysis."""

    # List of tuples: (Label, psi4.geometry string)
    molecule_defs = [
        ("Water", """
            0 1
            O
            H 1 0.957
            H 1 0.957 2 104.5
        """),
        ("Methanol", """
            0 1
            C
            H 1 1.09
            H 1 1.09 2 109.5
            H 1 1.09 2 109.5 3 120.0
            O 1 1.43 2 109.5 3 240.0
            H 5 0.96 4 109.5 1 60.0
        """),
        ("Carbon_Dioxide", """
            0 1
            C 
            O 1 1.16
            O 1 1.16 2 180.0
        """)
    ]
    return molecule_defs




# --- 4. MAIN EXECUTION FUNCTION (NEW) ---

def run_multiple_monomers_analysis():
    """Executes the analysis workflow for a list of small molecules."""

    # 1. Setup Environment
    setup_psi4_environment(METHOD, BASIS, THREADS)  # Increased memory for safety

    # 2. Define Geometries
    molecule_defs = define_molecules_for_test()

    all_results = {}

    # 3. Perform Calculations for each molecule
    for label, geom_string in molecule_defs:
        # Convert string to psi4.Molecule object
        print(geom_string)
        molecule = psi4.geometry(geom_string)

        # Run optimization and analysis
        E_opt, wfn_opt = optimize_and_get_energy(molecule, METHOD, BASIS, label)

        if E_opt != 999999.999999:
            calculate_properties_and_save(wfn_opt, label)
            all_results[label] = E_opt
        else:
            all_results[label] = "Failed"

    # 4. Final Summary
    my_log_print("\n" + "=" * 50)
    my_log_print("        ✅ Monomer Analysis Complete ✅")
    my_log_print("=" * 50)

    for label, energy in all_results.items():
        if isinstance(energy, float):
            my_log_print(f"Molecule: {label:<15} Energy: {energy:.8f} Hartree")
        else:
            my_log_print(f"Molecule: {label:<15} Status: {energy}")
    my_log_print("=" * 50)


if __name__ == '__main__':
    # Set the main log file name
    start = datetime.now()
    # We rename the main psi4 log file to match the new analysis
    psi4.set_output_file('monomer_analysis_detail.log', False)

    run_multiple_monomers_analysis()

    end = datetime.now()
    delta = end - start
    my_log_print(f"Total run time: {delta} seconds.")
    my_log_print(f"Threads count = {THREADS}")