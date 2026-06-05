"""
Rotational transform (iota) calculation and analysis module.

This module provides functions for:
- Calculating iota profile from fieldline data
- Detecting rational surfaces
- Adjusting LCFS to avoid rational surfaces
"""

import numpy as np
from scipy import stats


def calculate_iota_profile(X_lines, Y_lines, Z_lines, nturns, nphi, nlines):
    """
    Calculate the rotational transform iota profile from fieldline coordinates.
    
    Uses linear regression on poloidal angle θ vs toroidal angle φ to determine
    the rotational transform ι = dθ/dφ for each flux surface.
    
    Parameters
    ----------
    X_lines : ndarray
        Fieldline X coordinates, shape (nlines, nturns, nphi)
    Y_lines : ndarray
        Fieldline Y coordinates, shape (nlines, nturns, nphi)
    Z_lines : ndarray
        Fieldline Z coordinates, shape (nlines, nturns, nphi)
    nturns : int
        Number of toroidal turns
    nphi : int
        Number of poloidal points per turn
    nlines : int
        Total number of fieldlines
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'iota': array of rotational transform values
        - 'iota_err': array of fitting errors (standard error)
        - 'rho_mean': array of mean radial displacement from axis
        - 'rho_std': array of radial displacement standard deviation
    """
    # Calculate R coordinates
    R_lines = np.sqrt(X_lines**2 + Y_lines**2)
    
    # Magnetic axis coordinates
    axis_X = X_lines[0, :, :]
    axis_Y = Y_lines[0, :, :]
    axis_R = np.sqrt(axis_X**2 + axis_Y**2)
    axis_Z = Z_lines[0, :, :]
    
    # Calculate displacements from axis
    delt_r = R_lines[1:] - axis_R[np.newaxis, :, :]
    delt_r = np.reshape(delt_r, (nlines-1, nturns * nphi))
    delt_z = Z_lines[1:] - axis_Z[np.newaxis, :, :]
    delt_z = np.reshape(delt_z, (nlines-1, nturns * nphi))
    
    # Calculate poloidal angle
    theta_raw = np.arctan2(delt_z, delt_r)
    dtheta = np.diff(theta_raw, axis=1)
    dtheta[dtheta < -np.pi] += 2 * np.pi
    dtheta[dtheta > np.pi] -= 2 * np.pi
    theta_cumsum = np.cumsum(dtheta, axis=1)
    
    # Calculate radial displacement magnitude
    rho_raw = np.sqrt(delt_r**2 + delt_z**2)
    
    # Construct toroidal angle array
    phi_single = np.linspace(0, 2*np.pi, nphi, endpoint=False)
    PHI = np.tile(phi_single, nturns)
    PHI = np.tile(PHI, (nlines-1, 1))
    for i in range(nlines-1):
        PHI[i, :] += np.repeat(np.arange(nturns) * 2 * np.pi, nphi)
    
    # Pad theta to match PHI length
    theta_iota = np.hstack([np.zeros((nlines-1, 1)), theta_cumsum])
    
    # Initialize result arrays
    iota = np.zeros(nlines-1)
    iota_err = np.zeros(nlines-1)
    rho_mean = np.zeros(nlines-1)
    rho_std = np.zeros(nlines-1)
    
    # Perform linear fit for each fieldline
    for j in range(nlines-1):
        dex = R_lines[j+1, :, :].flatten() != 0
        
        if np.sum(dex) > 2:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                PHI[j, dex], theta_iota[j, dex]
            )
            
            iota[j] = slope
            iota_err[j] = std_err
            rho_mean[j] = np.mean(rho_raw[j, dex])
            rho_std[j] = np.std(rho_raw[j, dex])
        else:
            iota[j] = np.nan
            iota_err[j] = np.nan
            rho_mean[j] = np.nan
            rho_std[j] = np.nan
    
    return {
        'iota': iota,
        'iota_err': iota_err,
        'rho_mean': rho_mean,
        'rho_std': rho_std
    }


def check_rational_surface(iota, max_order=10, tolerance=0.01):
    """
    Check if iota value is close to a low-order rational surface.
    
    A rational surface occurs where ι = n/m, which can lead to island
    formation and reduced confinement.
    
    Parameters
    ----------
    iota : float
        Rotational transform value
    max_order : int, optional
        Maximum mode number (m, n) to check
    tolerance : float, optional
        Tolerance range for proximity detection
        
    Returns
    -------
    tuple or None
        If close to rational surface: (m, n, rational_value, error)
        Otherwise: None
    """
    closest_match = None
    min_error = float('inf')
    
    # Search all possible (m,n) combinations
    for m in range(1, max_order + 1):
        for n in range(1, max_order + 1):
            # Compute rational iota: ι = n/m
            rational_iota = n / m
            error = abs(rational_iota - iota)
            
            # If error < tolerance and best so far
            if error < tolerance and error < min_error:
                min_error = error
                closest_match = (m, n, rational_iota, error)
    
    return closest_match


def adjust_lcfs_avoid_rational_surface(lcfs_idx, iota, radius, 
                                      max_order=10, tolerance=0.01, 
                                      max_iterations=20):
    """
    Adjust LCFS index to avoid rational surfaces.
    
    Moves LCFS inward if it's too close to a low-order rational surface,
    which could cause problems for VMEC equilibrium calculations.
    
    Parameters
    ----------
    lcfs_idx : int
        Initial LCFS index
    iota : ndarray
        Rotational transform array
    radius : ndarray
        Radius array
    max_order : int, optional
        Maximum rational surface order to check
    tolerance : float, optional
        Tolerance for rational surface detection
    max_iterations : int, optional
        Maximum adjustment iterations
        
    Returns
    -------
    int
        Adjusted LCFS index
    """
    original_lcfs_idx = lcfs_idx
    iteration = 0
    
    print(f"Starting check for LCFS proximity to rational surfaces...")
    print(f"Initial index: {lcfs_idx}, radius: {radius[lcfs_idx]:.6f} m, iota: {iota[lcfs_idx]:.6f}")
    print("-" * 60)
    
    while lcfs_idx > 0 and iteration < max_iterations:
        check_iota = check_rational_surface(iota[lcfs_idx], max_order=max_order, tolerance=tolerance)
        
        if check_iota:
            m, n, rational_iota, error = check_iota
            print(f"Iteration {iteration+1}: LCFS (index={lcfs_idx}) close to rational surface ({m}, {n})")
            print(f"  Rational iota = {rational_iota:.6f}, current iota = {iota[lcfs_idx]:.6f}, error = {error:.6f}")
            lcfs_idx -= 1  # Move inward one layer
            iteration += 1
        else:
            # Found safe location
            print(f"✓ Found safe location (index={lcfs_idx})")
            break
    
    print("-" * 60)
    
    # Output adjustment results
    if iteration > 0:
        if lcfs_idx > 0:
            print(f"✓ LCFS adjusted: from index {original_lcfs_idx} -> {lcfs_idx} (moved inward {iteration} layers)")
            print(f"✓ New LCFS: radius={radius[lcfs_idx]:.6f} m, iota={iota[lcfs_idx]:.6f}")
        else:
            print(f"⚠ Warning: Reached innermost layer (index=0), cannot adjust further")
            lcfs_idx = 1  # At least preserve first layer
            print(f"✓ Forced LCFS index to {lcfs_idx}")
    else:
        print(f"✓ LCFS not close to rational surface, no adjustment needed (index={lcfs_idx})")
    
    print(f"\nFinal LCFS: index={lcfs_idx}, radius={radius[lcfs_idx]:.6f} m, iota={iota[lcfs_idx]:.6f}")
    
    return lcfs_idx

