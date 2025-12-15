from __future__ import print_function, absolute_import, division
from . import _fieldline_tracer
import f90wrap.runtime
import logging
import numpy

class Fieldline_Tracer(f90wrap.runtime.FortranModule):
    """
    Module fieldline_tracer
    
    
    Defined at src/fortran/fieldline_tracer_module.f90 lines 8-290
    
    """
    @staticmethod
    def set_verbose(verbose):
        """
        set_verbose(verbose)
        
        
        Defined at src/fortran/fieldline_tracer_module.f90 lines 45-48
        
        Parameters
        ----------
        verbose : bool
        
        """
        _fieldline_tracer.f90wrap_fieldline_tracer__set_verbose(verbose=verbose)
    
    @staticmethod
    def initialize_coils(coils_data):
        """
        initialize_coils(coils_data)
        
        
        Defined at src/fortran/fieldline_tracer_module.f90 lines 56-77
        
        Parameters
        ----------
        coils_data : float array
        
        """
        _fieldline_tracer.f90wrap_fieldline_tracer__initialize_coils(coils_data=coils_data)
    
    @staticmethod
    def cleanup_coils():
        """
        cleanup_coils()
        
        
        Defined at src/fortran/fieldline_tracer_module.f90 lines 83-92
        
        
        """
        _fieldline_tracer.f90wrap_fieldline_tracer__cleanup_coils()
    
    @staticmethod
    def trace_fieldlines(nturn, nphi, initial_rz, fieldline_data):
        """
        trace_fieldlines(nturn, nphi, initial_rz, fieldline_data)
        
        
        Defined at src/fortran/fieldline_tracer_module.f90 lines 104-161
        
        Parameters
        ----------
        nturn : int
        nphi : int
        initial_rz : float array
        fieldline_data : float array
        
        """
        _fieldline_tracer.f90wrap_fieldline_tracer__trace_fieldlines(nturn=nturn, \
            nphi=nphi, initial_rz=initial_rz, fieldline_data=fieldline_data)
    
    @staticmethod
    def get_magnetic_field(fieldpos, field):
        """
        get_magnetic_field(fieldpos, field)
        
        
        Defined at src/fortran/fieldline_tracer_module.f90 lines 283-290
        
        Parameters
        ----------
        fieldpos : float array
        field : float array
        
        """
        _fieldline_tracer.f90wrap_fieldline_tracer__get_magnetic_field(fieldpos=fieldpos, \
            field=field)
    
    _dt_array_initialisers = []
    

fieldline_tracer = Fieldline_Tracer()

