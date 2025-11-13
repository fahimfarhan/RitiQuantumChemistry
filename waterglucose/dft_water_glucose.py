import psi4
import numpy as np
import os
from datetime import datetime

# --- 1. CONFIGURATION ---
METHOD = 'B3LYP-D3'
BASIS = '6-31G(d)'
THREADS = os.cpu_count() or 4
LOG_FILE = 'glucose_dimer_summary.log'  # Custom file for summary prints


# --- 2. CUSTOM LOGGING FUNCTION ---
def my_log_print(message):
    """Prints the message to the console and appends it to the summary log file."""
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S] ")
    full_message = timestamp + message

    print(message)
    try:
        with open(LOG_FILE, 'a') as f:
            f.write(full_message + '\n')
    except Exception as e:
        print(f"Warning: Could not write to log file {LOG_FILE}: {e}")


# --- 3. SETUP FUNCTION ---
def setup_psi4_environment(method, basis, threads, memory='2 GB'):
    """Sets up Psi4 memory, threads, and global calculation options (increased memory)."""
    msg = f"Setting up Psi4: {method}/{basis} with {threads} threads. (27 atoms)"
    my_log_print(msg)

    psi4.set_memory(memory)
    psi4.set_num_threads(threads)

    psi4.set_options({
        'basis': basis,
        'd_convergence': 1e-7,  # Relax convergence slightly for large system
        'e_convergence': 1e-7,
        'guess': 'sad',
        'geom_maxiter': 10,  # Increased max steps
        'print': 2,
        'scf_type': 'df',
        'freeze_core': True,
    })


# --- 4. GEOMETRY DEFINITIONS ---

def define_dimer_geometry():
    """
    Defines the Water-Glucose complex geometry.
    Glucose coordinates (atoms 1-24) are followed by Water (atoms 25-27).
    The water is positioned near the O6-H bond.
    """
    my_log_print("Defining Water-Glucose complex geometry (27 atoms)...")

    # 0 1 defines the entire complex as neutral singlet
    # The whole block is one molecule for optimization
    dimer_geom = psi4.geometry("""
        0 1
        # --- Glucose (24 Atoms) ---
        O   0.94192   1.49258   0.31174  # O5
        C  -0.08830   0.72911  -0.34757  # C1
        H  -0.87108   1.23668  -0.89883
        C   0.45788  -0.65866  -0.87023  # C2
        O   0.62768  -1.50320   0.21855
        H   0.05737  -2.22744   0.05193
        H   1.30906  -0.56994  -1.39626
        C  -0.57503  -1.34110  -1.63666  # C3
        O  -1.62125  -0.42838  -1.87979
        H  -2.29653  -0.70631  -1.25878
        H  -0.29177  -2.04786  -2.42060
        C  -1.07720   0.01639  -1.07001  # C4
        O  -2.10091   0.67280  -1.80281
        H  -2.82525   0.17061  -1.42777
        H  -1.52042  -0.57963  -0.30182
        C  -0.50577   1.03366  -0.00762  # C5
        C  -1.54512   2.10309   0.50296  # C6
        O  -2.73033   1.61908   0.99963  # O6
        H  -3.26252   2.30810   0.74932  # H6
        H  -1.12788   2.76672   1.26620
        H  -1.80628   2.77977  -0.30189
        O   1.39611   0.37568  -1.14441  # O1
        H   1.53123   0.07604  -2.04423
        H   1.81053   0.10651   0.40237  # H1

        # --- Water (3 Atoms: O_W, H_W1, H_W2) ---
        # Water placed near the O6-H6 group for H-bonding
        O  -3.21733   0.88764   0.47953   # Atom 25
        H  -4.10271   1.21404   0.32354   # Atom 26 (forms H-bond with O6-H6)
        H  -2.86376   0.00741   0.37042   # Atom 27

        units angstrom
        symmetry c1
    """)
    return dimer_geom


def define_monomer_glucose_geometry():
    """Defines the isolated glucose monomer geometry."""
    # Extracted first 24 atoms from the dimer definition
    monomer_geom = psi4.geometry("""
        0 1
        O   0.94192   1.49258   0.31174
        C  -0.08830   0.72911  -0.34757
        H  -0.87108   1.23668  -0.89883
        C   0.45788  -0.65866  -0.87023
        O   0.62768  -1.50320   0.21855
        H   0.05737  -2.22744   0.05193
        H   1.30906  -0.56994  -1.39626
        C  -0.57503  -1.34110  -1.63666
        O  -1.62125  -0.42838  -1.87979
        H  -2.29653  -0.70631  -1.25878
        H  -0.29177  -2.04786  -2.42060
        C  -1.07720   0.01639  -1.07001
        O  -2.10091   0.67280  -1.80281
        H  -2.82525   0.17061  -1.42777
        H  -1.52042  -0.57963  -0.30182
        C  -0.50577   1.03366  -0.00762
        C  -1.54512   2.10309   0.50296
        O  -2.73033   1.61908   0.99963
        H  -3.26252   2.30810   0.74932
        H  -1.12788   2.76672   1.26620
        H  -1.80628   2.77977  -0.30189
        O   1.39611   0.37568  -1.14441
        H   1.53123   0.07604  -2.04423
        H   1.81053   0.10651   0.40237
        units angstrom
        symmetry c1
    """)
    return monomer_geom


# --- 5. CORE CALCULATION FUNCTIONS ---

def optimize_and_get_energy(molecule, method, basis, label):
    """Performs geometry optimization with error handling and returns the result."""
    my_log_print(f"\n🔬 Starting optimization for {label}...")

    E_opt = None
    wfn_opt = None

    try:
        # Run single-point for unoptimized Molden file (if dimer)
        if "Dimer" in label:
            my_log_print("Running initial single-point for UNOPTIMIZED structure output...")
            E_sp, wfn_sp = psi4.energy(f'{method}/{basis}', molecule=molecule, return_wfn=True)
            psi4.molden(wfn_sp, f'{label.lower()}_unoptimized.molden')

            # --- SCRATCH MANAGEMENT 1: Clean up after the SP energy calculation ---
            psi4.core.clean()
            my_log_print(f"Scratch cleaned after unoptimized SP energy.")

            my_log_print(f"Molden file for unoptimized {label} saved.")

        # Run the actual optimization
        E_opt, wfn_opt = psi4.optimize(f'{method}/{basis}', molecule=molecule, return_wfn=True)

        # --- SCRATCH MANAGEMENT 2: Clean up after the OPTIMIZATION ---
        psi4.core.clean()
        my_log_print(f"Scratch cleaned after {label} optimization.")

        # Psi4 returns None on failure (SCF or geometry convergence)
        if E_opt is None:
            raise RuntimeError("Psi4 optimization failed to converge.")

        my_log_print(f"✅ Optimized {label} Energy = {E_opt:.8f} Hartree")
        return E_opt, wfn_opt

    except Exception as e:
        my_log_print(f"🛑 FATAL CALCULATION ERROR for {label}: {e}")
        # Always attempt final cleanup even if calculation failed
        psi4.core.clean()
        return 999999.999999, None


def calculate_properties_and_save(wfn, label):
    """Calculates molecular properties and saves the optimized Molden file."""
    my_log_print(f"\nCalculating properties (Dipole and ESP Charges) for optimized {label}...")

    psi4.oeprop(wfn, 'DIPOLE', 'ESP_CHARGES')

    molden_file = f'{label.lower()}_optimized.molden'
    psi4.molden(wfn, molden_file)
    my_log_print(f"💾 Molden file saved as '{molden_file}'")


# --- 6. MAIN EXECUTION FUNCTION ---

def run_glucose_dimer_analysis():
    """Executes the full quantum chemistry workflow for the Water-Glucose complex."""

    # 0. Setup and Cleanup
    if os.path.exists(LOG_FILE):
        os.remove(LOG_FILE)
    setup_psi4_environment(METHOD, BASIS, THREADS)

    # 1. Define Geometries
    dimer_geom = define_dimer_geometry()
    glucose_geom = define_monomer_glucose_geometry()
    water_geom = psi4.geometry("""0 1\nO\nH 1 0.957\nH 1 0.957 2 104.5""")  # Simple water monomer

    # 2. Perform Calculations
    E_dimer_opt, wfn_dimer = optimize_and_get_energy(dimer_geom, METHOD, BASIS, "Water_Glucose_Dimer")

    # Monomer A (Glucose) Optimization - CRITICAL STEP
    E_glucose_opt, _ = optimize_and_get_energy(glucose_geom, METHOD, BASIS, "Glucose_Monomer")

    # Monomer B (Water) Optimization
    E_water_opt, _ = optimize_and_get_energy(water_geom, METHOD, BASIS, "Water_Monomer")

    # 3. Property Analysis and Output for Dimer
    calculate_properties_and_save(wfn_dimer, "Water_Glucose_Dimer")

    # 4. Calculate Binding Energy
    # E_binding = E(A...B) - E(A) - E(B)
    E_binding = E_dimer_opt - E_glucose_opt - E_water_opt
    E_binding_kcal = E_binding * psi4.constants.hartree2kcalmol

    # 5. Final Summary
    summary_message = f"""
        ================== Water-Glucose Final Summary ==================
        Method/Basis  : {METHOD} / {BASIS} (CPU Optimized, 27 Atoms)
        E(Optimized Dimer)  : {E_dimer_opt:.8f} Hartree
        E(Optimized Glucose): {E_glucose_opt:.8f} Hartree
        E(Optimized Water)  : {E_water_opt:.8f} Hartree

        🔹 Interaction Energy (Binding Energy) 🔹
        Interaction (Hartree) : {E_binding:.6f} Hartree
        Interaction (kcal/mol): {E_binding_kcal:.2f} kcal/mol
        =================================================================
    """
    my_log_print(summary_message)


if __name__ == '__main__':
    # Set the main Psi4 output file (detailed log)
    start = datetime.now()
    psi4.set_output_file('water_glucose_analysis.log', False)
    run_glucose_dimer_analysis()
    end = datetime.now()
    my_log_print(f"\nTotal time: {end - start} seconds.")
    pass