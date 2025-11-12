from datetime import datetime

import psi4
import numpy as np

# --- 1. CONFIGURATION ---
# Global variables for settings (or pass them via a dictionary/config object)
METHOD = 'B3LYP-D3'  # DFT with Grimme's dispersion correction
BASIS = '6-31G(d)'
THREADS = 4
LOG_FILE = 'analysis_summary.log' # Custom file for summary prints

def my_log_print(message):
    """Prints the message to the console and appends it to the summary log file."""
    # Add timestamp for better log tracking
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S] ")
    full_message = timestamp + message

    # 1. Print to console
    print(message)

    # 2. Append to log file
    try:
        with open(LOG_FILE, 'a') as f:
            f.write(full_message + '\n')
    except Exception as e:
        print(f"Warning: Could not write to log file {LOG_FILE}: {e}")

# --- 2. SETUP FUNCTION ---
def setup_psi4_environment(method, basis, threads, memory='1 GB'):
    """Sets up Psi4 memory, threads, and global calculation options."""
    my_log_print(f"Setting up Psi4: {method}/{basis} with {threads} threads.")
    psi4.set_memory(memory)
    psi4.set_num_threads(threads)

    psi4.set_options({
        'basis': basis,
        'd_convergence': 1e-8,
        'e_convergence': 1e-8,
        'guess': 'sad',
        'geom_maxiter': 50,
        'print': 2,
        'scf_type': 'df',  # CPU Optimization: Density Fitting
        'freeze_core': True,
    })


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


def optimize_and_get_energy(molecule, method, basis, label):
    """Performs geometry optimization and returns the optimized energy and wavefunction."""
    my_log_print(f"\n🔬 Starting optimization for {label}...")

    # Run a single-point energy calculation first to get the unoptimized wfn for Molden
    if "Dimer" in label:
        E_sp, wfn_sp = psi4.energy(f'{method}/{basis}', molecule=molecule, return_wfn=True)
        psi4.molden(wfn_sp, f'{label.lower()}_unoptimized.molden')
        my_log_print(f"Molden file for unoptimized {label} saved.")

    E_opt, wfn_opt = psi4.optimize(f'{method}/{basis}', molecule=molecule, return_wfn=True)
    my_log_print(f"✅ Optimized {label} Energy = {E_opt:.8f} Hartree")
    return E_opt, wfn_opt


def calculate_properties_and_save(wfn, label):
    """Calculates molecular properties and saves the optimized Molden file."""
    my_log_print(f"\nCalculating properties (Dipole and Charges) for optimized {label}...")

    # Calculate properties: Dipole Moment and ESP charges
    psi4.oeprop(wfn, 'DIPOLE', 'ESP_CHARGES')

    # Generate MOLDEN file for visualization
    molden_file = f'{label.lower()}_optimized.molden'
    psi4.molden(wfn, molden_file)
    my_log_print(f"💾 Molden file saved as '{molden_file}'")


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
    E_dimer_opt, wfn_dimer = optimize_and_get_energy(dimer_geom, METHOD, BASIS, "Dimer")
    E_monomer_opt, wfn_monomer = optimize_and_get_energy(monomer_geom, METHOD, BASIS, "Monomer")

    # 4. Property Analysis and Output for Dimer
    calculate_properties_and_save(wfn_dimer, "Dimer")

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


if __name__ == '__main__':
    # Set the main log file name
    psi4.set_output_file('water_dimer_analysis.log', False)
    run_water_dimer_analysis()