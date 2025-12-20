#!/usr/bin/env python3
"""
Example: Basic usage of coils2vmec for pwO device.

This is the refactored main execution script.
All function definitions are now in the coils2vmec package.
"""

import os
import sys
import numpy as np
from pathlib import Path

# Add package to path if not installed
script_dir = Path(__file__).parent
package_dir = script_dir.parent
sys.path.insert(0, str(package_dir))

from simsopt.field import BiotSavart, load_coils_from_makegrid_file
from simsopt.geo import SurfaceRZFourier, ToroidalFlux
from python import DescurConfig

# Import all needed functions from coils2vmec package
import python as c2v

# =============================================================================
# Configuration Parameters
# =============================================================================

# Device and file paths
device_name = "h1"
tag = "test"
coil_file = script_dir / 'coils.h1_measure_raw'

# Output directories
fieldline_dir = script_dir.parent / 'test' / 'fieldlines'
output_directory = fieldline_dir
vmec_output_dir = script_dir.parent / 'test' / 'vmec'

# Fieldline tracing parameters
nlines = 99
nphi = 360
nturn = 400
nfp = 1
mpol = 8
extcur = None
initial_rz = np.asfortranarray([1.32, 0])

# DESCUR parameters
config = DescurConfig()
nphi_descur = 60
config.mu = mpol
config.nv = nphi_descur
config.ftol = 2e-6
config.niter = 500

# Control flags
trace_flag = False      # Perform fieldline tracing
descur_flag = True     # Generate and run DESCUR
plot_flag = False      # Master switch for plotting

# =============================================================================
# Main Execution
# =============================================================================
print(f"Device: {device_name}")
print(f"VMEC output directory: {vmec_output_dir.name}")
print(f"Current run output directory: {output_directory.name}")
print(f"Coil file: {coil_file.name}")

# Load coils for toroidal flux calculation
coilpath = str(coil_file)
coils = load_coils_from_makegrid_file(filename=coilpath, order=20, ppp=36)

if extcur is not None:
    c2v.coils_with_extcur(coils, extcur)

bs_tf = BiotSavart(coils)

# Setup output paths
output_directory.mkdir(parents=True, exist_ok=True)
hdf5_file = output_directory / 'fieldlines_output.h5'

# =========================================================================
# Step 1: Fieldline Tracing
# =========================================================================

if trace_flag:
    print(f"\n{'='*60}")
    print("Step 1: Read coil file and trace fieldlines")
    print(f"{'='*60}")
    
    # Import fieldline tracer module
    try:
        from fieldline_tracer import fieldline_tracer
    except ImportError:
        print("Warning: Could not import fieldline_tracer module")
    
    # Read coil data
    coils_data = c2v.read_coils_file(str(coil_file), extcur=extcur, save_discrete=True)
    
    # Initialize coils in Fortran module
    fieldline_tracer.initialize_coils(coils_data)
    fieldline_tracer.set_verbose(False)
    
    print(f"\nStep 2: Parallel trace for {nlines} fieldlines")
    
    # Find magnetic axis
    c2v.find_axis(initial_rz, xtol=1e-10, max_iter=200)
    
    # Trace fieldlines in parallel
    fieldlines_data = c2v.trace_fieldlines_parallel(
        initial_guess=initial_rz,
        n_fieldlines=nlines,
        nturn=nturn,
        nphi=nphi,
        coils_data=coils_data,
        n_workers=None
    )
    
    # Report results
    print("\n✓ Fieldline tracing finished:")
    print(f"  Successful traces: {fieldlines_data['n_success']}/{fieldlines_data['n_total']}")
    print(f"  Magnetic axis: R={fieldlines_data['axis_rz'][0]:.6f}, Z={fieldlines_data['axis_rz'][1]:.6f}")
    print(f"  LCFS: R={fieldlines_data['lcfs_rz'][0]:.6f}, Z={fieldlines_data['lcfs_rz'][1]:.6f}")
    
    # Save to HDF5
    c2v.save_fieldlines_hdf5(fieldlines_data, str(hdf5_file), compress=True)

# =========================================================================
# Step 2: Data Processing and Iota Analysis
# =========================================================================

print("\u2713 Reshaping fieldline data")

# Load from HDF5
fieldlines_data = c2v.load_fieldlines_hdf5(str(hdf5_file))
nline = fieldlines_data['nline']
nphi = fieldlines_data['nphi']
nturn = fieldlines_data['nturn']

# Convert to numpy arrays
ALL_LINES_list = [fl for fl in fieldlines_data['fieldlines']]
LINES = np.array(ALL_LINES_list)
ALL_LINES = np.reshape(LINES, (nline, nturn, nphi, 4))

# Extract coordinates
X_lines = ALL_LINES[:, :, :, 0]
Y_lines = ALL_LINES[:, :, :, 1]
Z_lines = ALL_LINES[:, :, :, 2]
B_lines = ALL_LINES[:, :, :, 3]
R_lines = np.sqrt(X_lines**2 + Y_lines**2)
Phi_lines = np.arctan2(Y_lines, X_lines)
Phi_lines = np.mod(Phi_lines, 2*np.pi)

# Compute iota profile
iota_results = c2v.calculate_iota_profile(X_lines, Y_lines, Z_lines, nturn, nphi, nline)
iota = iota_results['iota']
iota_err = iota_results['iota_err']
rho_mean = iota_results['rho_mean']
rho_std = iota_results['rho_std']
radius = R_lines[1:, 0, 0]

# Find LCFS and islands
lcfs_result = c2v.find_lcfs_and_islands(
    radius, iota, 
    threshold_factor=15, 
    cooldown_factor=5.0, 
    plot_flag=False
)
lcfs_idx = lcfs_result['lcfs_index'] if lcfs_result is not None else None
lcfs_idx = lcfs_idx - 3
island_array = lcfs_result['island_array'] if lcfs_result is not None else []
ddiota_dr = lcfs_result['ddiota_dr'] if lcfs_result is not None else None

# Adjust LCFS to avoid rational surfaces
lcfs_idx = c2v.adjust_lcfs_avoid_rational_surface(
    lcfs_idx, iota, radius,
    max_order=10,
    tolerance=0.1,
    max_iterations=20
)

if lcfs_idx is not None:
    print(f"LCFS location: index={lcfs_idx}, radius={radius[lcfs_idx]:.6f} m, iota={iota[lcfs_idx]:.6f}")

# Optional plotting
if plot_flag:
    c2v.plot_iota_with_radius(iota, R_lines[1:, 0, 0], iota_err, lcfs_index=lcfs_idx)

# Fit iota profile
s = (np.linspace(0, rho_mean[lcfs_idx], lcfs_idx, endpoint=True) / rho_mean[lcfs_idx])**2
poly_order = 5
AI = np.polyfit(s, iota[:lcfs_idx], poly_order)[::-1]
print(f"Iota profile polynomial fit coefficients (order {poly_order}): {AI}")

Rlines_lcfs = R_lines[lcfs_idx]
print(f"LCFS fieldline R max: {Rlines_lcfs[0, 0]} m")

# More optional plotting
if plot_flag:
    fieldline_indices_to_plot = [0, lcfs_idx-1]
    c2v.plot_fieldlines_3d(X_lines, Y_lines, Z_lines, B_lines, 
                            fieldline_indices_to_plot, 
                            title='Selected Magnetic Fieldlines 3D')
    c2v.plot_poincare_sections(R_lines, Z_lines, phi_angles_deg=[0, 180], 
                                lcfs_index=lcfs_idx, POINT_SIZE=0.5, 
                                island_array=island_array, color_islands=False, 
                                title_suffix='(Colored by Flux Surface Region)')

# =========================================================================
# Step 3: DESCUR Surface Fitting
# =========================================================================

outcurve_path = output_directory / 'outcurve'

if descur_flag:
    print(f"✓ Generate and run DESCUR, output dir: {outcurve_path.name}")
    c2v.run_descur_python(
        R_lines=R_lines,
        Z_lines=Z_lines,
        Phi_lines=Phi_lines,
        lcfs_idx=lcfs_idx,
        nfp=1,
        nphi_descur=nphi_descur,
        output_directory=str(output_directory),
        config=config
    )
else:
    print(f"Skip DESCUR generation and run, check output dir: {outcurve_path.name}")

# Load fitted surface
surface = SurfaceRZFourier.from_vmec_input(str(outcurve_path))

# Optional self-intersection check
intersecting_check = False
if intersecting_check:
    print("Checking for self-intersection:")
    angles_rad = np.linspace(0, 2*np.pi, nphi_descur, endpoint=False)
    intersecting_angles = []
    for angle in angles_rad:
        is_intersecting_angle = surface.is_self_intersecting(angle=angle)
        if is_intersecting_angle:
            intersecting_angles.append(angle)
    if intersecting_angles:
        print(f"  Self-intersecting at angles: {[f'{a:.2f} rad ({a*180/np.pi:.1f}°)' for a in intersecting_angles]}")
    else:
        print("  No self-intersection detected at checked angles")

surface.get_quadpoints(180, 120)
if plot_flag:
    surface.plot(engine='plotly', plot_normal=True, close=True)

# =========================================================================
# Step 4: VMEC Input Generation
# =========================================================================

tf = ToroidalFlux(surface, bs_tf)
mpol = surface.mpol
ntor = surface.ntor
surface.stellsym = False
lasym = not surface.stellsym
phiedge = tf.J()
print(f"Toroidal flux at LCFS: {phiedge:.6e} Wb")

# Optional Poincare plot with surface
if plot_flag:
    phi_array = np.linspace(0, 360, 12, endpoint=False)
    output_path = Path(output_directory)
    c2v.plot_poincare_with_surface(
        surface, lcfs_idx, phi_array, R_lines, Z_lines, lcfs_idx, 
        figsize=(20, 60), plot_all=False
    )

# Save VMEC input file
c2v.save_vmec_input_for_surface(surface, bs_tf, extcur, mpol, ntor, lasym, nfp, vmec_output_dir, tag=tag, device_name=device_name, AI=AI, vmec_run=False)

print("\n✓ All steps completed successfully!")

