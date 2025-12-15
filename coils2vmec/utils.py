"""
Utility functions for VMEC and coil handling.

This module provides functions for:
- Setting coil currents
- Generating VMEC input files
- Post-processing DESCUR output
"""

import numpy as np
import os
from pathlib import Path
from simsopt.util import MpiPartition
from simsopt.geo import ToroidalFlux
from simsopt.mhd import Vmec


def coils_with_extcur(coils, extcur):
    """
    Set coil currents from extcur array.
    
    Parameters
    ----------
    coils : list
        List of coil objects
    extcur : array_like
        External current array
    """
    original_currents = [coil.current.get_value() for coil in coils]
    extcur_indices = [i for i in range(1, len(original_currents)) 
                    if original_currents[i] != original_currents[i-1]]
    
    if len(extcur) != len(extcur_indices) + 1:
        raise ValueError("extcur length does not match number of distinct current segments")

    j = 0
    extcur_current = extcur[j]
    for i in range(len(original_currents)):
        if i in extcur_indices:
            j = j + 1
            extcur_current = extcur[j] 
        coils[i].current.set_dofs([extcur_current])
    
    print("Coil currents set according to extcur.")


def save_vmec_input_for_surface(surface, bs_tf, extcur, mpol, ntor, lasym, nfp, 
                                save_path, tag, device_name="device", AI=None, 
                                vmec_run=False):
    """
    Save VMEC input file for fixed boundary equilibrium.
    
    Parameters
    ----------
    surface : simsopt surface object
        Boundary surface
    bs_tf : BiotSavart object
        Magnetic field object
    extcur : array_like
        External currents
    mpol : int
        Poloidal mode number
    ntor : int
        Toroidal mode number
    lasym : bool
        Whether to use asymmetric modes
    nfp : int
        Number of field periods
    save_path : Path
        Output directory
    tag : str
        Tag for filename
    device_name : str, optional
        Device name for filename
    AI : array_like, optional
        Iota profile coefficients
    vmec_run : bool, optional
        Whether to run VMEC after generating input
    """
    mpi = MpiPartition(ngroups=1)
    vmec = Vmec()
    vmec.indata.ns_array[:4] = [4, 9, 49, 99]
    vmec.indata.niter_array[:4] = [5000, 5000, 5000, 10000]
    vmec.indata.ftol_array[:4] = [1e-7, 1e-10, 1e-12, 1e-13]
    vmec.indata.pres_scale = 0
    
    if AI is not None:
        vmec.indata.ncurr = 0
        vmec.indata.ai[:len(AI)] = AI
        print(f"Enforce iota profile with AI: {AI}")
    else:
        vmec.indata.ncurr = 1
        vmec.indata.curtor = 0
        
    vmec.boundary = surface
    tf = ToroidalFlux(surface, bs_tf)
    tf_value = tf.J()
    
    if lasym:
        vmec.indata.lasym = lasym
    
    vmec.indata.extcur[:5] = extcur
    vmec.indata.phiedge = tf_value
    vmec.indata.mpol = mpol + 2
    vmec.indata.ntor = ntor + 2
    vmec.indata.ntheta = 4 * (mpol + 1)
    vmec.indata.nzeta = 4 * ntor
    vmec.indata.nfp = nfp
    
    # Fixed boundary
    vmec.indata.lfreeb = False
    fixed_name = f"input.{device_name}_{tag}_m{mpol}_n{ntor}_fixed"
    save_path = Path(save_path)
    save_path.mkdir(parents=True, exist_ok=True)
    vmec.write_input(save_path / fixed_name)
    print(f"Saved fixed boundary VMEC input file: {save_path / fixed_name}")
    
    if vmec_run:
        try:
            print("Running VMEC...")
            vmec.run()
            print("VMEC execution completed.")
        except Exception as e:
            print(f"VMEC execution failed: {str(e)}")


def ensure_indata_and_closure_in_outcurve(output_directory):
    """
    Ensure outcurve file has &indata before first RBC and '/' after last RBC.
    
    This is needed for VMEC compatibility when using Fortran DESCUR output.
    
    Parameters
    ----------
    output_directory : str
        Directory containing outcurve file
    """
    # Read original file
    outcurve_path = f"{output_directory}/outcurve"
    
    if not os.path.exists(outcurve_path):
        print(f"Warning: outcurve file not found at {outcurve_path}")
        return
    
    with open(outcurve_path, 'r') as file:
        lines = file.readlines()

    # Find first and last lines containing "RBC"
    rbc_idx = []
    modified = False

    for i, line in enumerate(lines):
        if ' RBC' in line:
            rbc_idx.append(i)
    
    first_rbc_idx = rbc_idx[1] if len(rbc_idx) > 1 else None
    last_rbc_idx = rbc_idx[-1] if rbc_idx else None

    # Add &indata before first RBC
    if first_rbc_idx is not None:
        if first_rbc_idx > 0 and '&indata' not in lines[first_rbc_idx-1]:
            lines.insert(first_rbc_idx, '&indata\n')
            print(f"✓ Added &indata before RBC at line {first_rbc_idx+1}")
            last_rbc_idx += 1
            modified = True
        elif first_rbc_idx == 0:
            lines.insert(0, '&indata\n')
            print(f"✓ RBC at file start, added &indata at beginning")
            last_rbc_idx += 1
            modified = True
        else:
            print(f"✓ &indata already present")

    # Add closure symbol after last RBC
    if last_rbc_idx is not None:
        if last_rbc_idx + 1 < len(lines):
            next_line = lines[last_rbc_idx + 1].strip()
            if not next_line.startswith('/'):
                lines.insert(last_rbc_idx + 1, '/\n')
                print(f"✓ Added closure '/' after last RBC at line {last_rbc_idx+1}")
                modified = True
            else:
                print(f"✓ Closure '/' already present")
        else:
            # Last RBC is at end of file
            lines.append('/\n')
            print(f"✓ Added closure '/' at end of file")
            modified = True

    # Write back file
    if modified:
        with open(outcurve_path, 'w') as f:
            f.writelines(lines)
        print(f"✓ Updated file: {outcurve_path}")


print("✓ Utility functions module loaded")
