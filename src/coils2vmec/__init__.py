"""
coils2vmec - Magnetic fieldline tracing and VMEC input generation toolkit.

This package provides tools for:
- Fieldline tracing from discrete coil configurations
- Rotational transform (iota) profile calculation
- LCFS and magnetic island detection
- DESCUR surface fitting
- VMEC input file generation
"""

__version__ = "0.1.0"

# Import main functions for easy access
from .fieldline import (
    initialize_coils,
    cleanup_coils,
    set_fortran_verbose,
    read_coils_file,
    find_axis,
    find_lcfs,
    trace_fieldlines_parallel,
    save_fieldlines_hdf5,
    load_fieldlines_hdf5
)

from .iota import (
    calculate_iota_profile,
    check_rational_surface,
    adjust_lcfs_avoid_rational_surface
)

from .lcfs import (
    find_lcfs_and_islands
)

from .descur import (
    save_lcfs_for_descur,
    run_descur_python,
    run_descur_silent,
    generate_and_run_descur
)

from .descur_python import (
    DescurConfig,
    DescurFitter
)

from .utils import (
    coils_with_extcur,
    save_vmec_input_for_surface,
    ensure_indata_and_closure_in_outcurve
)

from .plotting import (
    plot_axis_3d,
    plot_iota_with_radius,
    plot_fieldlines_3d,
    plot_poincare_sections,
    plot_surface_cross_section_RZ,
    plot_surface_cross_sections_multi,
    plot_poincare_with_surface
)

__all__ = [
    # Fieldline functions
    'initialize_coils',
    'cleanup_coils',
    'set_fortran_verbose',
    'read_coils_file',
    'find_axis',
    'find_lcfs',
    'trace_fieldlines_parallel',
    'save_fieldlines_hdf5',
    'load_fieldlines_hdf5',
    
    # Iota analysis
    'calculate_iota_profile',
    'check_rational_surface',
    'adjust_lcfs_avoid_rational_surface',
    
    # LCFS detection
    'find_lcfs_and_islands',
    
    # DESCUR
    'save_lcfs_for_descur',
    'run_descur_python',
    'run_descur_silent',
    'generate_and_run_descur',
    'DescurConfig',
    'DescurFitter',
    
    # VMEC utilities
    'coils_with_extcur',
    'save_vmec_input_for_surface',
    'ensure_indata_and_closure_in_outcurve',
    
    # Plotting
    'plot_axis_3d',
    'plot_iota_with_radius',
    'plot_fieldlines_3d',
    'plot_poincare_sections',
    'plot_surface_cross_section_RZ',
    'plot_surface_cross_sections_multi',
    'plot_poincare_with_surface',
]
