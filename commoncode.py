import os
from datetime import datetime

import psi4
import psutil # <-- New Import
from psi4.core import Molecule

# --- 1. CONFIGURATION ---
# Global variables for settings (or pass them via a dictionary/config object)
METHOD = 'B3LYP-D3'  # DFT with Grimme's dispersion correction
BASIS = '6-31G(d)'
THREADS = os.cpu_count() or 4

def my_log_print(message):
    """Prints the message to the console and appends it to the summary log file."""
    # Add timestamp for better log tracking
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S] ")
    full_message = timestamp + message

    # 1. Print to console
    print(message)
    LOG_FILE = 'analysis_summary.log'  # Custom file for summary prints
    # 2. Append to log file
    try:
        with open(LOG_FILE, 'a') as f:
            f.write(full_message + '\n')
    except Exception as e:
        print(f"Warning: Could not write to log file {LOG_FILE}: {e}")


def get_optimal_memory_setting(target_fraction=0.75, minimum_gb=4):
    """
    Reads total system RAM using psutil, calculates a safe fraction,
    and returns the Psi4 memory string (e.g., '24 GB').
    """
    try:
        # psutil.virtual_memory().total returns RAM in bytes
        total_ram_bytes = psutil.virtual_memory().total
        total_ram_gib = total_ram_bytes / (1024 ** 3)  # Convert to GiB

        # Calculate the target memory
        target_ram_gib = total_ram_gib * target_fraction

        # Ensure the memory is at least the minimum required for a meaningful run
        final_ram_gib = max(target_ram_gib, minimum_gb)

        # Convert to an integer value in GB (which is often what clusters expect)
        # We use a slightly safer integer Gigabyte conversion for Psi4's string input
        memory_string = f"{int(final_ram_gib)} GB"

        my_log_print(f"Detected System RAM: {total_ram_gib:.2f} GiB.")
        my_log_print(f"Setting Psi4 memory to: {memory_string} (approx. {target_fraction * 100:.0f}% of total).")
        return memory_string

    except Exception as e:
        # Fallback if psutil fails (e.g., on a restricted HPC environment)
        my_log_print(f"Warning: Could not read system RAM ({e}). Falling back to minimum setting.")
        return f"{minimum_gb} GB"

def setup_psi4_environment(method, basis, threads, memory=None):
    """Sets up Psi4 memory, threads, and global calculation options."""
    # Use the dynamic memory setting if no memory is provided (None)
    if memory is None:
        memory = get_optimal_memory_setting()

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
    pass


def optimize_and_get_energy(molecule: Molecule, method, basis, label):
    """Performs geometry optimization with error handling, returns result,
       and handles complex-specific scratch/molden output via fragment count."""

    my_log_print(f"\n🔬 Starting optimization for {label}...")

    E_opt = None
    wfn_opt = None

    # Check if the molecule is a complex (i.e., defined with multiple fragments via '--')
    is_complex = len(molecule.get_fragments()) > 1

    try:
        # --- COMPLEX-SPECIFIC LOGIC (For Dimers, Trimers, etc.) ---
        if is_complex:
            my_log_print("Running initial single-point for UNOPTIMIZED complex structure output...")

            # 1. Single-Point calculation
            E_sp, wfn_sp = psi4.energy(f'{method}/{basis}', molecule=molecule, return_wfn=True)

            # 2. Output Unoptimized Molden
            psi4.molden(wfn_sp, f'{label.lower()}_unoptimized.molden')
            my_log_print(f"Molden file for unoptimized complex saved.")

            # 3. SCRATCH MANAGEMENT 1: Clean up after the SP energy calculation
            psi4.core.clean()
            my_log_print(f"Scratch cleaned after unoptimized SP energy.")

            # Note: The molecule object already holds the initial geometry for the optimization below

        # --- UNIVERSAL OPTIMIZATION STEP ---
        my_log_print(f"Starting geometry optimization for {label}...")
        E_opt, wfn_opt = psi4.optimize(f'{method}/{basis}', molecule=molecule, return_wfn=True)

        # --- SCRATCH MANAGEMENT 2: Clean up after the OPTIMIZATION ---
        psi4.core.clean()
        my_log_print(f"Scratch cleaned after {label} optimization.")

        # Psi4 returns None on failure (SCF or geometry convergence)
        if E_opt is None:
            # Check the detailed output for the reason (e.g., failed to converge)
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
    if wfn is None:
        my_log_print(f"\nThe wave function for {label} is Null. cannot compute optimization! returning...")
        return

    my_log_print(f"\nCalculating properties (Dipole and Charges) for optimized {label}...")

    # Calculate properties: Dipole Moment and ESP charges
    psi4.oeprop(wfn, 'DIPOLE', 'ESP_CHARGES')

    # Generate MOLDEN file for visualization
    molden_file = f'{label.lower()}_optimized.molden'
    psi4.molden(wfn, molden_file)
    my_log_print(f"💾 Molden file saved as '{molden_file}'")
    psi4.core.clean()
    pass
