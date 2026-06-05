"""
LCFS (Last Closed Flux Surface) detection module.

This module provides functions for detecting the LCFS and magnetic islands
using second derivative analysis of the iota profile.
"""

import numpy as np
import matplotlib.pyplot as plt


def find_lcfs_and_islands(radius, iota, threshold_factor=5.0, 
                          cooldown_factor=5.0, plot_flag=False):
    """
    Detect LCFS and island locations by searching for d²iota/dr² peaks.
    
    The algorithm searches for peaks in the second derivative of iota with
    respect to radius. The first high peak is identified as LCFS using a
    stricter threshold (3x island threshold). Subsequent peaks indicate
    island boundaries. A cooldown mechanism prevents double-counting nearby peaks.
    
    Parameters
    ----------
    radius : ndarray
        Radial coordinate array
    iota : ndarray
        Rotational transform array
    threshold_factor : float, optional
        Island threshold multiplier (relative to global median).
        LCFS uses 3x this factor.
    cooldown_factor : float, optional
        Cooldown period as percentage of array length to avoid
        detecting multiple peaks in the same structure
    plot_flag : bool, optional
        Whether to produce diagnostic plots
        
    Returns
    -------
    dict
        Contains:
        - 'lcfs_index': LCFS index
        - 'lcfs_radius': LCFS radius value
        - 'island_array': List of island boundary indices
        - 'median_value': Global median of |d²iota/dr²|
        - 'threshold_island': Island detection threshold
        - 'threshold_lcfs': LCFS detection threshold
        - 'ddiota_dr': Second derivative array
        - 'cooldown_length': Cooldown period in array indices
    """
        
    # Check input array lengths
    if len(radius) != len(iota) or len(radius) < 2:
        print("Error: Input arrays must have the same length and at least 2 elements.")
        return {
            'lcfs_index': -1, 'lcfs_radius': np.nan, 'island_array': [], 
            'median_value': np.nan, 'threshold_island': np.nan, 
            'threshold_lcfs': np.nan, 'ddiota_dr': None
        }
        
    # Calculate second derivative d(diota/dr)/dr
    diota_dr = np.gradient(iota, radius)
    ddiota_dr = np.gradient(diota_dr, radius)
    abs_ddiota = np.abs(ddiota_dr)

    # Establish global robust baseline and thresholds
    mid_abs_ddiota = np.median(abs_ddiota)
    threshold_island = mid_abs_ddiota * threshold_factor  # island threshold (lenient)
    threshold_lcfs = mid_abs_ddiota * threshold_factor * 3.0  # LCFS threshold (strict, 3x island)
    
    # Calculate cooldown length (cooldown_factor as percentage)
    cooldown_length = max(1, int(len(radius) * cooldown_factor / 100.0))
    
    # Initialize result storage and state flags
    island_array = []
    lcfs_index = None
    # State flag: True means currently in a detected "peak" region
    is_in_peak = False
    # Cooldown counter: tracks distance from last peak detection
    cooldown_counter = 0
    
    # Search for peak region start points from axis outward
    for i in range(len(radius)):
        current_value = abs_ddiota[i]
        
        if current_value > threshold_island:
            # Current point exceeds threshold
            if lcfs_index is None and current_value > threshold_lcfs:
                lcfs_index = i-1  # First peak is the LCFS

            if not is_in_peak and cooldown_counter == 0:
                # Current point not in peak region and not in cooldown
                # This marks start of a new peak region
                island_array.append(i)
                is_in_peak = True  # Set flag, mark entry to peak region
                
        else:
            # Current point below threshold
            if is_in_peak:
                # Exiting a peak region, reset flag and start cooldown
                is_in_peak = False
                cooldown_counter = cooldown_length  # Start cooldown
        
        # Update cooldown counter
        if cooldown_counter > 0:
            cooldown_counter -= 1

    # Determine LCFS and finalize results
    if lcfs_index is None:
        # No peaks found, use last point as default LCFS
        lcfs_index = len(radius) - 1
        
    lcfs_radius = radius[lcfs_index]
        
    lcfs_result = {
        'lcfs_index': lcfs_index,
        'lcfs_radius': lcfs_radius,
        'island_array': island_array,
        'median_value': mid_abs_ddiota,
        'threshold_island': threshold_island,
        'threshold_lcfs': threshold_lcfs,
        'ddiota_dr': ddiota_dr,
        'cooldown_length': cooldown_length
    }
    
    # Plotting section (controlled by plot_flag)
    if plot_flag:
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        # --- Primary Axis: LCFS Metric ---
        color_ddiota = 'tab:blue'
        ax1.set_xlabel('Radius (r)')
        ax1.set_ylabel(r'$|d^2 \iota / dr^2|$ (LCFS Metric)', color=color_ddiota)
        
        # Plot |d²iota/dr²|
        ax1.plot(radius, abs_ddiota, label=r'$|d^2 \iota / dr^2|$', color=color_ddiota)
        ax1.tick_params(axis='y', labelcolor=color_ddiota)
        
        # Plot global median and thresholds
        ax1.axhline(mid_abs_ddiota, color='green', linestyle='--', 
                    label=f'Global Median: {mid_abs_ddiota:.3e}')
        ax1.axhline(threshold_island, color='orange', linestyle='-.', linewidth=1.5,
                    label=f'Island Threshold ({threshold_factor}x Median): {threshold_island:.3e}')
        ax1.axhline(threshold_lcfs, color='red', linestyle='-.', 
                    label=f'LCFS Threshold ({threshold_factor*3:.1f}x Median): {threshold_lcfs:.3e}')
        
        # --- Secondary Axis: Iota Profile ---
        ax2 = ax1.twinx() 
        color_iota = 'tab:orange'
        ax2.set_ylabel(r'$\iota$ (Rotational Transform)', color=color_iota)
        ax2.plot(radius, iota, label=r'$\iota$ Profile', color=color_iota, linestyle=':')
        ax2.tick_params(axis='y', labelcolor=color_iota)

        # --- Mark LCFS and Island Boundaries ---
        
        # Mark LCFS
        if island_array:
            lcfs_r = lcfs_result['lcfs_radius']
            lcfs_idx = lcfs_result['lcfs_index']
            
            ax1.axvline(lcfs_r, color='purple', linestyle='--', linewidth=2, 
                        label=f'LCFS Start at r = {lcfs_r:.3f}')
            ax1.scatter(lcfs_r, abs_ddiota[lcfs_idx], 
                        color='purple', marker='o', s=100, zorder=5, label='LCFS Start Point')
            ax2.scatter(lcfs_r, iota[lcfs_idx], 
                        color='purple', marker='x', s=100, zorder=5)

        # Mark other Island/Structure boundaries
        if len(island_array) > 1:
            for i, idx in enumerate(island_array[1:]):
                r_i = radius[idx]
                ax1.axvline(r_i, color='gray', linestyle=':', linewidth=1)
                ax1.scatter(r_i, abs_ddiota[idx], color='gray', marker='^', s=50, zorder=4, 
                            label=f'Island {i+1}' if i == 0 else None)

        # Organize legend and title
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='best', fontsize='small')
        
        plt.title(f'LCFS and Island Detection (LCFS: {threshold_factor*3:.1f}x, Island: {threshold_factor}x Median, Cooldown: {cooldown_factor}% = {cooldown_length} pts)')
        ax1.grid(True)
        fig.tight_layout() 
        plt.show()

    return lcfs_result

