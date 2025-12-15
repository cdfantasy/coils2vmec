"""
Plotting module for magnetic field visualization.

This module provides functions for:
- Poincaré section plots
- 3D fieldline visualization
- Iota profile plots
- Surface cross-section plots
"""

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from pathlib import Path


def plot_axis_3d(axis_X, axis_Y, axis_Z):
    """
    Plot the magnetic axis as a 3D Plotly figure.
    
    Parameters
    ----------
    axis_X : ndarray
        X coordinates of magnetic axis
    axis_Y : ndarray
        Y coordinates of magnetic axis
    axis_Z : ndarray
        Z coordinates of magnetic axis
    """
    fig = go.Figure(data=[go.Scatter3d(
        x=axis_X.flatten(),
        y=axis_Y.flatten(),
        z=axis_Z.flatten(),
        mode='lines',
        line=dict(width=2, color='blue'),
        name='Magnetic Axis'
    )])

    fig.update_layout(
        title='Magnetic Axis Fieldline',
        scene=dict(
            xaxis_title='X (m)',
            yaxis_title='Y (m)',
            zaxis_title='Z (m)',
            aspectmode='data'
        ),
        width=800,
        height=600
    )
    fig.show()


def plot_iota_with_radius(iota, Radius, iota_err=None, lcfs_index=None):
    """
    Plot iota profile as a function of radius with optional LCFS marker.
    
    Parameters
    ----------
    iota : ndarray
        Rotational transform values
    Radius : ndarray
        Radius values
    iota_err : ndarray, optional
        Iota fitting errors
    lcfs_index : int, optional
        LCFS index for marking
    """
    fig, ax1 = plt.subplots(figsize=(8, 6))
    
    # Primary Y-axis: iota
    color1 = 'blue'
    ax1.set_xlabel('Radius (m)')
    ax1.set_ylabel('iota', color=color1)
    ax1.scatter(Radius, iota, s=10, c=color1, alpha=0.7, label='iota')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, alpha=0.3)
    
    # Mark LCFS if provided
    if lcfs_index is not None:
        lcfs_r = Radius[lcfs_index]
        lcfs_iota = iota[lcfs_index]
        ax1.scatter([lcfs_r], [lcfs_iota], color='red', s=20, zorder=5, 
                   marker='o', facecolors='none', edgecolors='red', linewidths=2.5,
                   label=f'LCFS at r={lcfs_r:.4f}m')
    
    # Secondary Y-axis for error if provided
    if iota_err is not None:
        ax2 = ax1.twinx()
        color2 = 'red'
        ax2.set_ylabel('iota error', color=color2)
        ax2.plot(Radius, iota_err, color=color2, alpha=0.5, linewidth=2, label='iota error')
        ax2.tick_params(axis='y', labelcolor=color2)
    
    ax1.legend(loc='upper left')
    plt.title('Rotational Transform Profile')
    fig.tight_layout()
    plt.show()


def plot_fieldlines_3d(X_lines, Y_lines, Z_lines, B_lines, fieldline_indices, 
                       title='Magnetic Fieldlines 3D (B-Field Colored)'):
    """
    Plot selected fieldlines in 3D colored by magnetic field strength.
    
    Parameters
    ----------
    X_lines : ndarray
        X coordinates of fieldlines
    Y_lines : ndarray
        Y coordinates of fieldlines
    Z_lines : ndarray
        Z coordinates of fieldlines
    B_lines : ndarray
        Magnetic field strength |B| for coloring
    fieldline_indices : list of int
        Indices of fieldlines to plot
    title : str, optional
        Plot title
    """
    fig = go.Figure()

    for idx in fieldline_indices:
        # Flatten coordinate and field data
        x_data = X_lines[idx].flatten()
        y_data = Y_lines[idx].flatten()
        z_data = Z_lines[idx].flatten()
        b_data = B_lines[idx].flatten()
        
        # Ensure data lengths match
        if not (len(x_data) == len(y_data) == len(z_data) == len(b_data)):
            print(f"Warning: Data lengths mismatch for fieldline {idx}. Skipping.")
            continue

        # Use Scatter3d with field-based coloring
        fig.add_trace(go.Scatter3d(
            x=x_data,
            y=y_data,
            z=z_data,
            mode='lines',
            line=dict(
                width=5,
                color=b_data,
                colorscale='Viridis',
                colorbar=dict(title='|B| (T)')
            ),
            name=f'Fieldline {idx+1}',
            showlegend=False 
        ))

    fig.update_layout(
        title={
            'text': title,
            'y':0.9,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        scene=dict(
            xaxis_title='X (m)',
            yaxis_title='Y (m)',
            zaxis_title='Z (m)',
            aspectmode='data' 
        ),
        width=800,
        height=600
    )

    fig.show()


def plot_poincare_sections(R_lines, Z_lines, phi_angles_deg, lcfs_index, POINT_SIZE, 
                          island_array=None, color_islands=True, title_suffix=''):
    """
    Plot Poincaré sections (R, Z) at specified toroidal angles, colored by flux surface region.

    Parameters
    ----------
    R_lines : ndarray
        Major radius coordinates (nlines, nturns, nphi)
    Z_lines : ndarray
        Z coordinates (nlines, nturns, nphi)
    phi_angles_deg : list/ndarray
        Toroidal angles in degrees for Poincaré sections
    lcfs_index : int
        Index of LCFS
    POINT_SIZE : float
        Marker size for Poincaré points
    island_array : list of int, optional
        Indices of island boundaries
    color_islands : bool, optional
        Whether to color distinct island regions
    title_suffix : str, optional
        Suffix for plot title
    """
    # Total fieldlines (axis + surfaces)
    nlines_total = R_lines.shape[0] 
    # Surfaces only
    nlines_surfaces = nlines_total - 1 
    
    # Convert input angles to radians
    phi_angles_rad = np.deg2rad(phi_angles_deg)
    
    # --- 1. Determine color classification ---
    
    # Assign color group to each non-axis fieldline
    color_groups = np.zeros(nlines_surfaces, dtype=int)
    
    # Color groups: 0=Axis, 1=Core, 2=LCFS, 3=SOL/Stochastic
    
    if 0 <= lcfs_index < nlines_surfaces:
        color_groups[:lcfs_index] = 1 # Core
        color_groups[lcfs_index] = 2 # LCFS
        color_groups[lcfs_index + 1:] = 3 # SOL
    else:
        color_groups[:] = 1 

    # --- Island/Structure coloring (if enabled) ---
    if color_islands and island_array and len(island_array) > 0:
        structures = [lcfs_index] + [idx for idx in island_array if idx > lcfs_index]
        structures = sorted(list(set(structures)))
        
        current_color_idx = 4 # Start island region coloring at 4
        
        for i in range(len(structures)):
            start_idx = structures[i] 
            end_idx = structures[i+1] if i + 1 < len(structures) else nlines_surfaces
            
            if start_idx >= lcfs_index:
                 color_groups[start_idx:end_idx] = current_color_idx
            
            if start_idx == lcfs_index:
                 color_groups[lcfs_index] = 2 

            if i == len(structures) - 1 and end_idx < nlines_surfaces:
                 color_groups[end_idx:] = current_color_idx 
            
            current_color_idx += 1
            
    # Define base color map (Plotly colors)
    colors_base = ['darkblue', 'black', 'red', 'gray']
    
    if color_islands and np.max(color_groups) >= len(colors_base):
        max_group = np.max(color_groups)
        # Use Plotly HSV colorscale for additional colors
        additional_colors = px.colors.sample_colorscale("HSV", max_group - len(colors_base) + 2)
        colors_base.extend(additional_colors[1:])
            
    # --- 2. Setup plotting ---
    
    nphi_data = R_lines.shape[2]
    phi_data_rad = np.linspace(0, 2*np.pi, nphi_data, endpoint=False)
    
    num_plots = len(phi_angles_rad)
    cols = min(2, num_plots)
    rows = int(np.ceil(num_plots / cols))

    # Create subplots
    fig = make_subplots(
        rows=rows, cols=cols,
        subplot_titles=[f'φ = {deg:.1f}°' for deg in phi_angles_deg],
        vertical_spacing=0.12,
        horizontal_spacing=0.10
    )

    # --- 3. Plot Poincaré sections ---
    
    for plot_idx, phi_target_rad in enumerate(phi_angles_rad):
        phi_target_deg = phi_angles_deg[plot_idx]
        row = plot_idx // cols + 1
        col = plot_idx % cols + 1

        idx_in_phi = np.argmin(np.abs(phi_data_rad - (phi_target_rad % (2 * np.pi))))
        
        legend_labels = {
            0: 'Magnetic Axis',
            1: 'Core',
            2: 'LCFS',
            3: 'SOL / Stochastic'
        }
        
        # Plot magnetic axis
        fig.add_trace(
            go.Scatter(
                x=R_lines[0, :, idx_in_phi].flatten(),
                y=Z_lines[0, :, idx_in_phi].flatten(),
                mode='markers',
                marker=dict(
                    symbol='x',
                    size=4,
                    color=colors_base[0]
                ),
                name=legend_labels[0],
                legendgroup='axis',
                showlegend=(plot_idx == 0)
            ),
            row=row, col=col
        )
        
        # Plot all flux surfaces
        added_groups = set()
        
        for i in range(nlines_surfaces):
            group_id = color_groups[i]
            line_color = colors_base[group_id]
            
            # Determine legend label
            if group_id > 3:
                label_name = f'Structure/Island Region {group_id-3}'
            else:
                label_name = legend_labels.get(group_id, '')
            
            # Show in legend only first time encountering group
            show_in_legend = (plot_idx == 0 and group_id not in added_groups and label_name != '')
            
            # Plot R vs Z
            fig.add_trace(
                go.Scatter(
                    x=R_lines[i+1, :, idx_in_phi].flatten(),
                    y=Z_lines[i+1, :, idx_in_phi].flatten(),
                    mode='markers',
                    marker=dict(
                        size=POINT_SIZE,
                        color=line_color
                    ),
                    name=label_name,
                    legendgroup=f'group_{group_id}',
                    showlegend=show_in_legend
                ),
                row=row, col=col
            )
            
            if show_in_legend:
                added_groups.add(group_id)

        # Set axis properties
        fig.update_xaxes(
            title_text='R (m)',
            scaleanchor=f'y{plot_idx+1}',
            scaleratio=1,
            row=row, col=col
        )
        fig.update_yaxes(
            title_text='Z (m)',
            row=row, col=col
        )

    # Set overall layout
    fig.update_layout(
        title={
            'text': f'Poincaré Sections Colored by Flux Surface Region {title_suffix}',
            'y':0.98,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        height=400 * rows,
        width=500 * cols,
        showlegend=True,
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        )
    )
    
    fig.show()


def plot_surface_cross_section_RZ(surface, phi=0, ntheta=100, ax=None, **plot_kwargs):
    """
    Plot surface cross-section in R-Z plane at given phi angle.
    
    Parameters
    ----------
    surface : simsopt surface object
        Surface to plot
    phi : float, optional
        Poloidal angle (0-1, where 0=0°, 0.5=180°, 1=360°)
    ntheta : int, optional
        Number of poloidal points
    ax : matplotlib axis, optional
        Axis to plot on
    **plot_kwargs : optional
        Arguments for plt.plot()
    
    Returns
    -------
    tuple
        (R, Z) coordinates of cross-section
    """
    # Get cross-section points
    cross_pts = surface.cross_section(phi, thetas=ntheta)
    R = np.sqrt(cross_pts[:, 0]**2 + cross_pts[:, 1]**2)
    Z = cross_pts[:, 2]
    
    if ax is None:
        fig, ax = plt.subplots()
    
    ax.plot(R, Z, **plot_kwargs)
    ax.set_xlabel('R')
    ax.set_ylabel('Z')
    ax.set_title(f'Cross-section at phi={phi*360:.1f}°')
    ax.axis('equal')
    plt.show()
    
    return R, Z


def plot_surface_cross_sections_multi(surface, angles_deg=[0, 60, 120, 180, 240, 300], 
                                      ntheta=100, figsize=(15, 10)):
    """
    Plot surface cross-sections at multiple angles.
    
    Parameters
    ----------
    surface : simsopt surface object
        Surface to plot
    angles_deg : list, optional
        List of angles in degrees
    ntheta : int, optional
        Number of poloidal points
    figsize : tuple, optional
        Figure size
        
    Returns
    -------
    tuple
        (fig, axes) matplotlib objects
    """
    # Create subplots
    fig_clumn = len(angles_deg) // 2 + len(angles_deg) % 2
    fig, axes = plt.subplots(2, fig_clumn, figsize=figsize)
    axes = axes.flatten()
    
    for i, angle_deg in enumerate(angles_deg):
        # Convert angle from degrees to phi (0-1)
        phi = (angle_deg / 360.0) % 1.0
        
        # Get cross-section points
        cross_pts = surface.cross_section(phi, thetas=ntheta)
        R = np.sqrt(cross_pts[:, 0]**2 + cross_pts[:, 1]**2)
        Z = cross_pts[:, 2]
        
        # Plot to corresponding subplot
        axes[i].plot(R, Z, 'b-', linewidth=2)
        axes[i].set_xlabel('R (m)')
        axes[i].set_ylabel('Z (m)')
        axes[i].set_title(f'φ = {angle_deg}° (φ={phi:.2f})')
        axes[i].axis('equal')
        axes[i].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return fig, axes


def plot_poincare_with_surface(surface, lcfs_idx, phi_array, R_matrix, Z_matrix, lcfs, 
                               figsize=(20, 20), plot_all=True, tf=None, mpol=None, 
                               ntor=None, output_path=None):
    """
    Display Poincaré plots at multiple phi angles with surface overlay.
    
    Parameters
    ----------
    surface : simsopt surface object
        Fitted surface to overlay
    lcfs_idx : int
        LCFS index
    phi_array : array_like
        Toroidal angles in degrees
    R_matrix : ndarray
        R coordinates (nlines, nphi, nturns)
    Z_matrix : ndarray
        Z coordinates (nlines, nphi, nturns)
    lcfs : int
        LCFS line index
    figsize : tuple, optional
        Figure size
    plot_all : bool, optional
        Whether to plot all points or just LCFS
    tf : ToroidalFlux object, optional
        For title
    mpol : int, optional
        For title and filename
    ntor : int, optional
        For title and filename
    output_path : Path, optional
        Save path for figure
        
    Returns
    -------
    None
    """
    nlines, nphi, nturns = R_matrix.shape
    n_plots = len(phi_array)
    
    # Calculate subplot layout
    cols = min(2, n_plots)
    rows = (n_plots + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    if n_plots == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    if tf is not None and mpol is not None and ntor is not None:
        fig.suptitle(f'Surface and Poincaré at toroidal flux = {tf.J():.3e} with m = {mpol}, n = {ntor}', 
                    fontsize=16)

    for i, phi_deg in enumerate(phi_array):
        ax = axes[i]
        
        # Convert angle to index
        phi_deg = int(phi_deg) % 360
        phi_idx = phi_deg / 180 * np.pi
        cross_pts = surface.cross_section(phi_idx, thetas=100)
        R = np.sqrt(cross_pts[:, 0]**2 + cross_pts[:, 1]**2)
        Z = cross_pts[:, 2]
        ax.plot(R, Z, 'b-', linewidth=2, label='Surface')
        
        if not plot_all:
            closest_R = R_matrix[lcfs_idx, :, phi_deg]
            closest_Z = Z_matrix[lcfs_idx, :, phi_deg]
            ax.scatter(closest_R, closest_Z, s=1, color='black', marker='o', label='Closest Point')
        else:
            # Extract all points at this phi angle
            R_phi = R_matrix[:, :, phi_deg]
            Z_phi = Z_matrix[:, :, phi_deg]
            # Flatten data
            R_flat = R_phi.flatten()
            Z_flat = Z_phi.flatten()
            
            # Plot all fieldline points
            if i == 0:
                print(f'plot all points, from {R_phi[0,0]:.3e} m to {R_phi[-1,0]:.3e} m')
            ax.scatter(R_flat, Z_flat, s=0.5, alpha=0.5, color='black')
        
        # Mark LCFS
        if lcfs < nlines:
            R_lcfs = R_matrix[lcfs, :, phi_deg]
            Z_lcfs = Z_matrix[lcfs, :, phi_deg]
            ax.scatter(R_lcfs, Z_lcfs, s=2, color='red', marker='x', label='LCFS')
            ax.legend()
        
        ax.set_xlabel('R (m)')
        ax.set_ylabel('Z (m)')
        ax.set_title(f'φ = {phi_deg}°')
        ax.axis('equal')
        ax.grid(True, alpha=0.3)
    
    # Hide extra subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    
    if output_path is not None and mpol is not None and ntor is not None:
        plt_savepath = Path(output_path)
        plt.savefig(plt_savepath / f'm{mpol}_n{ntor}.png', dpi=300)
    
    plt.show()


print("✓ Plotting module loaded")
