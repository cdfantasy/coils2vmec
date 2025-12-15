"""
DESCUR (surface fitting) module.

This module provides functions for fitting magnetic flux surfaces using
DESCUR (Discrete Equilibrium Surface CURve fitting) algorithm.
"""

import numpy as np
import os
import subprocess


def save_lcfs_for_descur(R_lines, Z_lines, Phi_lines, lcfs_idx, nfp=1, 
                         nphi_descur=72, output_directory='.'):
    """
    Save LCFS data in DESCUR format.
    
    Extracts LCFS data from fieldline arrays and formats for DESCUR input.
    Uses greedy nearest-neighbor algorithm to reorder points for smooth curve.
    
    Parameters
    ----------
    R_lines : ndarray
        R coordinate array, shape (nlines, nturns, nphi)
    Z_lines : ndarray
        Z coordinate array, shape (nlines, nturns, nphi)
    Phi_lines : ndarray
        Phi coordinate array, shape (nlines, nturns, nphi)
    lcfs_idx : int
        LCFS fieldline index
    nfp : int, optional
        Number of field periods
    nphi_descur : int, optional
        Number of toroidal angle points for DESCUR
    output_directory : str, optional
        Output directory path
        
    Returns
    -------
    None
        Writes 'descur_input' file to output_directory
    """
    # Build complete file path
    output_path = os.path.join(output_directory, 'descur_input')
    
    # Extract LCFS data and transpose
    R_lcfs = R_lines[lcfs_idx].T
    Z_lcfs = Z_lines[lcfs_idx].T
    phi_lcfs = Phi_lines[lcfs_idx].T
    
    R_slice = R_lcfs[0, :]
    Z_slice = Z_lcfs[0, :]

    nturns = R_lcfs.shape[1]
    # Initialize sort indices and visited flags
    sorted_indices = np.zeros(nturns, dtype=int)
    visited = np.zeros(nturns, dtype=bool)

    # Start from first point
    current_index = 0
    sorted_indices[0] = current_index
    visited[current_index] = True

    # Greedy nearest-neighbor search
    for k in range(1, nturns):
        current_R = R_slice[current_index]
        current_Z = Z_slice[current_index]
        
        min_dist_sq = np.inf
        next_index = -1
        
        # Find nearest unvisited point
        for i in range(nturns):
            if not visited[i]:
                dist_sq = (R_slice[i] - current_R)**2 + (Z_slice[i] - current_Z)**2
                
                if dist_sq < min_dist_sq:
                    min_dist_sq = dist_sq
                    next_index = i
        
        # Update current point
        if next_index != -1:
            current_index = next_index
            sorted_indices[k] = current_index
            visited[current_index] = True
        else:
            break

    # Apply reordering
    R_lcfs_reordered = R_lcfs[:, sorted_indices]
    Z_lcfs_reordered = Z_lcfs[:, sorted_indices]
    phi_lcfs_reordered = phi_lcfs[:, sorted_indices]

    # Ensure first toroidal angle is 0
    phi_lcfs_reordered[0, :] = 0
    
    # Get total toroidal angle points
    nphi = R_lcfs.shape[0]
    
    # Select toroidal angle indices for output
    descur_idx = np.linspace(0, nphi, nphi_descur, endpoint=False).astype(int)
    
    # Extract corresponding data
    R_lcfs_descur = R_lcfs_reordered[descur_idx]
    Z_lcfs_descur = Z_lcfs_reordered[descur_idx]
    phi_lcfs_descur = phi_lcfs_reordered[descur_idx]
    
    # Prepare DESCUR data format: R, phi, Z
    descur_data = np.column_stack((R_lcfs_descur.flatten(), 
                                  phi_lcfs_descur.flatten(), 
                                  Z_lcfs_descur.flatten()))
    
    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)
    
    # Write file
    with open(output_path, 'w') as f:
        f.write(f'{nturns} {nphi_descur} {nfp}\n')
        np.savetxt(f, descur_data, fmt='%.10e')
    
    print(f"DESCUR input file saved: {output_path}")


def run_descur_silent(working_directory='.'):
    """
    Run xdescur silently with hardcoded filename 'descur_input'.
    
    Parameters
    ----------
    working_directory : str, optional
        Working directory containing descur_input file
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    # Switch to working directory
    original_dir = os.getcwd()
    os.chdir(working_directory)
    
    commands = [
        "4",      # spectral convergence parameter
        "v",      # VMEC-compatible compression
        "0",      # read data from file
        "descur_input"  # hardcoded filename
    ]
    
    input_text = "\n".join(commands) + "\n"
    
    try:
        process = subprocess.Popen(
            ['xdescur'],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True
        )
        
        process.communicate(input=input_text)
        
        # Switch back to original directory
        os.chdir(original_dir)
        return True
        
    except FileNotFoundError:
        print("Error: xdescur command not found")
        os.chdir(original_dir)
        return False
    except Exception as e:
        print(f"Error running xdescur: {e}")
        os.chdir(original_dir)
        return False


def generate_and_run_descur(R_lines, Z_lines, Phi_lines, lcfs_idx, 
                           nfp=1, nphi_descur=72, 
                           output_directory='.', run=True):
    """
    Generate DESCUR input file and execute xdescur silently.
    
    Parameters
    ----------
    R_lines : ndarray
        R coordinate array
    Z_lines : ndarray
        Z coordinate array
    Phi_lines : ndarray
        Phi coordinate array
    lcfs_idx : int
        LCFS index
    nfp : int, optional
        Number of field periods
    nphi_descur : int, optional
        Number of toroidal points
    output_directory : str, optional
        Output directory path
    run : bool, optional
        Whether to run xdescur after generating input
        
    Returns
    -------
    bool or None
        If run=True, returns success status; otherwise None
    """
    # 1. Generate DESCUR input file
    save_lcfs_for_descur(
        R_lines=R_lines,
        Z_lines=Z_lines, 
        Phi_lines=Phi_lines,
        lcfs_idx=lcfs_idx,
        nfp=nfp,
        nphi_descur=nphi_descur,
        output_directory=output_directory
    )
    
    # 2. Run xdescur silently
    if run:
        return run_descur_silent(working_directory=output_directory)


def run_descur_python(R_lines, Z_lines, Phi_lines, lcfs_idx, 
                     nfp=1, nphi_descur=72, 
                     output_directory='.', 
                     use_python=True,
                     log_filename='descur.log',
                     config=None):
    """
    Process fieldline data with Python DESCUR directly, no file I/O.
    
    Parameters
    ----------
    R_lines : ndarray
        R coordinate array, shape (nlines, nturns, nphi)
    Z_lines : ndarray
        Z coordinate array, shape (nlines, nturns, nphi)
    Phi_lines : ndarray
        Phi coordinate array, shape (nlines, nturns, nphi)
    lcfs_idx : int
        LCFS fieldline index
    nfp : int, optional
        Number of field periods
    nphi_descur : int, optional
        Number of toroidal angle points
    output_directory : str, optional
        Output directory path
    use_python : bool, optional
        Whether to use Python implementation (default True)
    log_filename : str, optional
        Log file name
    config : DescurConfig, optional
        DESCUR configuration object
        
    Returns
    -------
    dict or None
        Fitting results dictionary (if use_python=True), otherwise None
    """
    if use_python:
        # Use Python DESCUR implementation
        from .descur_python import DescurFitter, DescurConfig
        
        # Create default config if not provided
        if config is None:
            config = DescurConfig(
                mu=14,
                nv=nphi_descur,
                ftol=1e-5,
                niter=500,
                mexp=4,
                pexp=4.0,
                qexp=1.0
            )
        else:
            # If config provided, ensure nv matches nphi_descur
            config.nv = nphi_descur
        
        # Create fitter
        fitter = DescurFitter(config)
        
        # Set log file path
        log_path = os.path.join(output_directory, log_filename)
        fitter.setup_logger(log_file=log_path, console=True)
        
        # Prepare input from fieldline data (no file I/O)
        fitter.logger.info(f"\nProcessing fieldline data with Python DESCUR...")
        rin, zin = fitter.prepare_from_fieldlines(
            R_lines, Z_lines, Phi_lines, 
            lcfs_idx=lcfs_idx,
            nfp=nfp,
            nphi_descur=nphi_descur
        )
        
        # Execute fit
        fitter.logger.info("\nStarting DESCUR fit...")
        results = fitter.fit(rin, zin, log_file=log_path)
        
        # Write output file
        os.makedirs(output_directory, exist_ok=True)
        output_file = os.path.join(output_directory, 'outcurve')
        fitter.write_output(results, output_file)
        
        fitter.logger.info(f"\n✓ Python DESCUR finished, results saved to: {output_file}")
        fitter.logger.info(f"✓ Log saved to: {log_path}")
        
        return results
    else:
        # Use Fortran DESCUR (requires file I/O)
        print(f"\nUsing Fortran DESCUR (xdescur)...")
        generate_and_run_descur(
            R_lines, Z_lines, Phi_lines, lcfs_idx,
            nfp=nfp, nphi_descur=nphi_descur,
            output_directory=output_directory,
            run=True
        )
        
        # Post-process outcurve file
        from .utils import ensure_indata_and_closure_in_outcurve
        ensure_indata_and_closure_in_outcurve(output_directory)
        return None


print("✓ DESCUR module loaded")
