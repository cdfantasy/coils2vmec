"""
Fieldline tracing module for magnetic field analysis.

This module provides functions for:
- Reading coil files
- Finding magnetic axis
- Finding LCFS
- Parallel fieldline tracing
- Saving/loading fieldline data in HDF5 format
"""

import numpy as np
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import h5py
from scipy.optimize import root

# Import fieldline_tracer module from package
try:
    from .fieldline_tracer import Fieldline_Tracer
    Fieldline_Tracer.set_verbose(False)
except ImportError:
    Fieldline_Tracer = None
    print("Warning: Fieldline_Tracer module not found. Fieldline tracing functions will not work.")

def initialize_coils(coils_data):
    """
    Initialize coils in the Fortran module.
    
    Parameters
    ----------
    coils_data : ndarray
        Discrete coil data array (n_points x 4): [x, y, z, current]
    """
    if Fieldline_Tracer is not None:
        Fieldline_Tracer.initialize_coils(coils_data)
    else:
        raise ImportError("Fieldline_Tracer module is not available. Cannot initialize coils.")

def cleanup_coils():
    """
    Clean up coils in the Fortran module.
    """
    if Fieldline_Tracer is not None:
        Fieldline_Tracer.cleanup_coils()
    else:
        raise ImportError("Fieldline_Tracer module is not available. Cannot clean up coils.")

def set_fortran_verbose(verbose):
    """
    Control verbose output from the Fortran module.
    
    Parameters
    ----------
    verbose : bool
        Whether to enable verbose output
    """
    if Fieldline_Tracer is not None:
        Fieldline_Tracer.set_verbose(verbose)


def read_coils_file(filename, extcur=None, save_discrete=False):
    """
    Read coil file in coils.* format and convert to discrete coil format.
    
    Parameters
    ----------
    filename : str
        Path to coil file
    extcur : dict, optional
        External current dictionary {group_id: current_value}
        If provided, use specified currents instead of file values
    save_discrete : bool, optional
        Whether to save converted discrete_coil.txt file
        
    Returns
    -------
    ndarray
        Discrete coil data array (n_points x 4): [x, y, z, current]
    """
    print(f"\nReading coils file: {filename}")
    
    coils_group_idx = [0]
    line_num = 0
    group_id = 1
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith('#'):
                continue
            if line.startswith('periods') or line.startswith('begin') or line.startswith('mirror'):
                continue
            if line.startswith('end'):
                break

            line_num += 1
            try:
                parts = line.split()
                if len(parts) >= 5:
                    if group_id == int(parts[4]):
                        coils_group_idx[group_id-1] = line_num
                    else:
                        group_id = int(parts[4])
                        coils_group_idx.append(line_num)
            except (ValueError, IndexError):
                    pass
        
        print(f"  Coil File Analysis:")
        print(f"    Total coil groups: {len(coils_group_idx)}")
        print(f"    Coil group line indices: {coils_group_idx}")

    line_num = 0
    data_points = []
    extcur_id = 0
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith('#'):
                continue
            if line.startswith('periods') or line.startswith('begin') or line.startswith('mirror'):
                continue
            if line.startswith('end'):
                break

            line_num += 1
            parts = line.split()
            if len(parts) == 4:
                x = float(parts[0])
                y = float(parts[1])
                z = float(parts[2])
                current = float(parts[3])

                for i in range(len(coils_group_idx)):
                    if line_num >= coils_group_idx[i]:
                        extcur_id = i+1
                # Determine actual current to use
                if extcur is not None and extcur_id in extcur:
                    actual_current = extcur[extcur_id]
                else:
                    actual_current = current
                data_points.append([x, y, z, actual_current])

            elif len(parts) >= 5:
                x = float(parts[0])
                y = float(parts[1])
                z = float(parts[2])
                current = float(parts[3])
                data_points.append([x, y, z, current])

    discrete_coil = np.array(data_points, dtype=np.float64)
    
    print(f"  Total data points read: {len(discrete_coil)}")
    print(f"  Data range:")
    print(f"    x: [{discrete_coil[:,0].min():.4f}, {discrete_coil[:,0].max():.4f}]")
    print(f"    y: [{discrete_coil[:,1].min():.4f}, {discrete_coil[:,1].max():.4f}]")
    print(f"    z: [{discrete_coil[:,2].min():.4f}, {discrete_coil[:,2].max():.4f}]")
    print(f"    current: [{discrete_coil[:,3].min():.2e}, {discrete_coil[:,3].max():.2e}]")
    
    # Save discrete coil file
    if save_discrete:
        output_file = 'discrete_coil.txt'
        np.savetxt(output_file, discrete_coil, fmt='%.8e')
        print(f"  Saved discrete coil data to: {output_file}")
    
    return discrete_coil


def find_axis(initial_guess, xtol=1e-10, max_iter=200):
    """
    Find the magnetic axis by ensuring the fieldline returns to the start after one toroidal turn.
    
    Parameters
    ----------
    initial_guess : array_like
        Initial (R, Z) guess for the magnetic axis
    xtol : float, optional
        Tolerance for convergence
    max_iter : int, optional
        Maximum number of iterations
        
    Returns
    -------
    ndarray
        Magnetic axis position [R, Z]
    """
    print(f"\nSearching for magnetic axis...")
    print(f"  Initial guess: R={initial_guess[0]:.6f}, Z={initial_guess[1]:.6f}")
    
    def axis_residual(rz):
        n_points = 360
        initial_rz_array = np.array(rz, dtype=np.float64, order='F')
        fieldline_data = np.zeros((2*n_points, 4), dtype=np.float64, order='F')
        
        try:
            Fieldline_Tracer.trace_fieldlines(2, n_points, initial_rz_array, fieldline_data)
            final_xyz = fieldline_data[n_points, :3]
            final_R = np.sqrt(final_xyz[0]**2 + final_xyz[1]**2)
            final_Z = final_xyz[2]
            residual = np.array([final_R - rz[0], final_Z - rz[1]])
            return residual
        except Exception as e:
            print(f"  Warning: Fieldline tracing failed at R={rz[0]:.6f}, Z={rz[1]:.6f}: {e}")
            return np.array([1e10, 1e10])
    
    result = root(
        axis_residual,
        initial_guess,
        method='hybr',
        tol=xtol,
        options={'maxfev': max_iter * (len(initial_guess) + 1)}
    )
    
    final_residual = result.fun
    distance = np.linalg.norm(final_residual)
    
    print("  Optimization completed:")
    print(f"    Axis position: R={result.x[0]:.10f}, Z={result.x[1]:.10f}")
    print(f"    Distance error: {distance:.2e}")
    print(f"    Converged: {result.success}")
    
    return result.x


def find_lcfs(initial_guess, precision_order=1e-3, nturn=40, verbose=True):
    """
    Find the LCFS by outward binary search in R while keeping Z fixed.
    
    Parameters
    ----------
    initial_guess : array_like
        Initial (R, Z) guess for LCFS search
    precision_order : float, optional
        Target precision for binary search
    nturn : int, optional
        Number of toroidal turns for escape detection
    verbose : bool, optional
        Whether to print progress information
        
    Returns
    -------
    ndarray
        LCFS position [R, Z]
    """
    if verbose:
        print(f"\nSearching for LCFS...")
        print(f"  Initial guess: R={initial_guess[0]:.6f}, Z={initial_guess[1]:.6f}")
        print(f"  Target precision: {precision_order}")
    
    fixed_Z = initial_guess[1]
    initial_R = initial_guess[0]
    
    def check_escape(R):
        nphi = 30
        initial_rz = np.array([R, fixed_Z], dtype=np.float64, order='F')
        fieldline_data = np.zeros((nturn*nphi, 4), dtype=np.float64, order='F')
        
        try:
            Fieldline_Tracer.trace_fieldlines(nturn, nphi, initial_rz, fieldline_data)
            last_points = fieldline_data[-5:, :]
            if np.all(np.abs(last_points) < 1e-10):
                return True
            else:
                return False
        except Exception:
            return True
    
    if verbose:
        print("\n  Step 1: Finding escape boundary...")
    
    R_lower = initial_R
    R_upper = None
    escaped = check_escape(R_lower)
    
    if escaped:
        if verbose:
            print("  Error: Initial position already escaped!")
        return np.array([initial_R, fixed_Z])
    
    if verbose:
        print(f"  Initial position: R={R_lower:.6f}, tracing SUCCESS")
    
    search_R = R_lower
    search_step = 0.05
    max_search_R = 3.0
    search_iter = 0
    
    while search_R < max_search_R and search_iter < 50:
        search_R += search_step
        escaped = check_escape(search_R)
        search_iter += 1
        
        if escaped:
            R_upper = search_R
            if verbose:
                print(f"  Found escape at R={R_upper:.6f}")
            break
        else:
            if verbose and search_iter % 5 == 0:
                print(f"    Step {search_iter}: R={search_R:.6f}, tracing SUCCESS")
    
    if R_upper is None:
        R_upper = max_search_R
    
    if verbose:
        print(f"  Search bounds: R_lower={R_lower:.6f}, R_upper={R_upper:.6f}")
        print("\n  Step 2: Binary search for LCFS...")
    
    iterations = 0
    max_iter = int(np.log2((R_upper - R_lower) / precision_order)) + 10
    
    while (R_upper - R_lower) > precision_order and iterations < max_iter:
        iterations += 1
        R_mid = (R_lower + R_upper) / 2.0
        escaped = check_escape(R_mid)
        
        if escaped:
            R_upper = R_mid
        else:
            R_lower = R_mid
    
    lcfs_R = R_lower - 2*precision_order
    lcfs_rz = np.array([lcfs_R, fixed_Z])
    
    if verbose:
        print("\n  LCFS found:")
        print(f"    Position: R={lcfs_R:.10f}, Z={fixed_Z:.10f}")
        print(f"    Iterations: {iterations}")
    
    return lcfs_rz


def trace_single_fieldline(args):
    """
    Helper to trace a single fieldline (for parallelization).
    
    Parameters
    ----------
    args : tuple
        (index, rz, nturn, nphi, coils_data)
        
    Returns
    -------
    tuple
        (index, fieldline_data, success, error_message)
    """
    index, rz, nturn, nphi, coils_data = args
    
    set_fortran_verbose(False)
    Fieldline_Tracer.initialize_coils(coils_data)
    
    n_points = nturn * nphi
    initial_rz = np.array(rz, dtype=np.float64, order='F')
    fieldline_data = np.zeros((n_points, 4), dtype=np.float64, order='F')
    
    try:
        Fieldline_Tracer.trace_fieldlines(nturn, nphi, initial_rz, fieldline_data)
        return (index, fieldline_data, True, None)
    except Exception as e:
        return (index, None, False, str(e))
    finally:
        Fieldline_Tracer.cleanup_coils()


def trace_fieldlines_parallel(axis_rz, lcfs_rz, n_fieldlines, nturn=100, nphi=360, 
                              coils_data=None, n_workers=None):
    """
    Uniformly sample flux surfaces between axis and LCFS and trace in parallel.
    
    Parameters
    ----------
    axis_rz : array_like
        (R, Z) position of the magnetic axis
    lcfs_rz : array_like
        (R, Z) position of the LCFS
    n_fieldlines : int
        Number of fieldlines to trace
    nturn : int, optional
        Number of toroidal turns per fieldline
    nphi : int, optional
        Number of poloidal points per turn
    coils_data : ndarray
        Coil data array from read_coils_file
    n_workers : int, optional
        Number of parallel workers (None = auto)
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'fieldlines': list of fieldline data arrays
        - 'initial_positions': starting (R, Z) for each fieldline
        - 'indices': successful fieldline indices
        - 'n_total': total number of fieldlines requested
        - 'n_success': number of successfully traced fieldlines
        - 'execution_time': total execution time
        - 'axis_rz': magnetic axis position
        - 'lcfs_rz': LCFS position
        - 'nline', 'nphi', 'nturn': tracing parameters
    """
    print(f"\nTracing {n_fieldlines} fieldlines in parallel...")
    print(f"  Axis: R={axis_rz[0]:.6f}, Z={axis_rz[1]:.6f}")
    print(f"  LCFS: R={lcfs_rz[0]:.6f}, Z={lcfs_rz[1]:.6f}")
    
    if coils_data is None:
        raise ValueError("coils_data must be provided")
    
    direction = lcfs_rz - axis_rz
    lcfs_radius = np.linalg.norm(direction)
    direction = direction / lcfs_radius
    
    radii = np.linspace(0, lcfs_radius, n_fieldlines)
    initial_positions = [axis_rz + direction * r for r in radii]
    
    print(f"\nTracing {len(initial_positions)} fieldlines (nturn={nturn}, nphi={nphi})...")
    
    tasks = [(i, rz, nturn, nphi, coils_data) for i, rz in enumerate(initial_positions)]
    
    t_start = time.time()
    results = [None] * n_fieldlines
    success_flags = [False] * n_fieldlines
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(trace_single_fieldline, task): task[0] for task in tasks}
        
        for future in as_completed(futures):
            index, data, success, error = future.result()
            results[index] = data
            success_flags[index] = success
            
            if (index + 1) % max(1, n_fieldlines // 10) == 0:
                print(f"  Progress: {sum(success_flags)}/{n_fieldlines} completed")
    
    t_end = time.time()
    elapsed = t_end - t_start
    
    # Remove failed fieldlines
    successful_fieldlines = []
    successful_positions = []
    successful_indices = []
    
    for i in range(n_fieldlines):
        if success_flags[i]:
            successful_fieldlines.append(results[i])
            successful_positions.append(initial_positions[i])
            successful_indices.append(i)
    
    n_success = len(successful_fieldlines)
    total_points = n_success * nturn * nphi
    nline = n_success

    print("\nFieldline tracing completed:")
    print(f"  Execution time: {elapsed:.4f} seconds")
    print(f"  Successful traces: {n_success}/{n_fieldlines}")
    print(f"  Total points traced: {total_points}")
    if n_success > 0:
        print(f"  Throughput: {total_points/elapsed:.1f} points/second")
    
    return {
        'fieldlines': successful_fieldlines,
        'initial_positions': successful_positions,
        'indices': successful_indices,
        'n_total': n_fieldlines,
        'n_success': n_success,
        'execution_time': elapsed,
        'axis_rz': axis_rz,
        'lcfs_rz': lcfs_rz,
        'nline': n_fieldlines,
        'nphi': nphi,
        'nturn': nturn
    }


def save_fieldlines_hdf5(fieldlines_data, filename, compress=True):
    """
    Save fieldline data to HDF5 format.
    
    Parameters
    ----------
    fieldlines_data : dict
        Fieldline data dictionary from trace_fieldlines_parallel
    filename : str
        Output HDF5 filename
    compress : bool, optional
        Whether to use gzip compression
    """
    print(f"\nSaving fieldlines to HDF5: {filename}")
    
    with h5py.File(filename, 'w') as f:
        fieldlines_group = f.create_group('fieldlines')
        for i, fieldline in enumerate(fieldlines_data['fieldlines']):
            if compress:
                fieldlines_group.create_dataset(
                    f'fieldline_{i:04d}', 
                    data=fieldline,
                    compression='gzip',
                    compression_opts=4
                )
            else:
                fieldlines_group.create_dataset(f'fieldline_{i:04d}', data=fieldline)
        
        initial_pos_array = np.array(fieldlines_data['initial_positions'])
        f.create_dataset('initial_positions', data=initial_pos_array)
        f.create_dataset('indices', data=np.array(fieldlines_data['indices']))
        
        f.attrs['n_total'] = fieldlines_data['n_total']
        f.attrs['n_success'] = fieldlines_data['n_success']
        f.attrs['axis_R'] = fieldlines_data['axis_rz'][0]
        f.attrs['axis_Z'] = fieldlines_data['axis_rz'][1]
        f.attrs['lcfs_R'] = fieldlines_data['lcfs_rz'][0]
        f.attrs['lcfs_Z'] = fieldlines_data['lcfs_rz'][1]
        f.attrs['nline'] = fieldlines_data['nline']
        f.attrs['nphi'] = fieldlines_data['nphi']
        f.attrs['nturn'] = fieldlines_data['nturn']
        
    
    print(f"  ✓ Saved {fieldlines_data['n_success']} fieldlines to {filename}")
    if os.path.exists(filename):
        print(f"  File size: {os.path.getsize(filename) / 1024 / 1024:.2f} MB")


def load_fieldlines_hdf5(filename):
    """
    Load fieldline data from HDF5 file.
    
    Parameters
    ----------
    filename : str
        Input HDF5 filename
        
    Returns
    -------
    dict
        Fieldline data dictionary
    """
    print(f"\nLoading fieldlines from HDF5: {filename}")
    
    with h5py.File(filename, 'r') as f:
        fieldlines = []
        fieldlines_group = f['fieldlines']
        n_fieldlines = len(fieldlines_group.keys())
        
        for i in range(n_fieldlines):
            fieldline = fieldlines_group[f'fieldline_{i:04d}'][:]
            fieldlines.append(fieldline)
        
        initial_positions = f['initial_positions'][:].tolist()
        indices = f['indices'][:].tolist()
        
        n_total = f.attrs['n_total']
        n_success = f.attrs['n_success']
        nphi = f.attrs['nphi']
        nturn = f.attrs['nturn']
        nline = f.attrs['nline']
        axis_rz = np.array([f.attrs['axis_R'], f.attrs['axis_Z']])
        lcfs_rz = np.array([f.attrs['lcfs_R'], f.attrs['lcfs_Z']])
    
    print(f"  ✓ Loaded {n_success} fieldlines from {filename}")
    
    return {
        'fieldlines': fieldlines,
        'initial_positions': initial_positions,
        'indices': indices,
        'nline': nline,
        'nphi': nphi,
        'nturn': nturn,
        'axis_rz': axis_rz,
        'lcfs_rz': lcfs_rz
    }

