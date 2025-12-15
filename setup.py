#!/usr/bin/env python3
"""
Setup script for coils2vmec package with Fortran compilation
"""

import os
import subprocess
import sys
from pathlib import Path
from setuptools import setup
from setuptools.command.build_ext import build_ext

class FortranBuildExt(build_ext):
    """Custom build_ext command to compile Fortran using Makefile"""

    def run(self):
        # First, build Fortran modules using Makefile
        print("\n" + "="*70)
        print("Building Fortran extensions using Makefile...")
        print("="*70 + "\n")
        
        try:
            # Check for required tools
            subprocess.run(['gfortran', '--version'], 
                         capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("ERROR: gfortran not found. Please install a Fortran compiler.")
            print("  Ubuntu/Debian: sudo apt-get install gfortran")
            print("  macOS: brew install gcc")
            print("  CentOS/RHEL: sudo yum install gcc-gfortran")
            sys.exit(1)

        try:
            subprocess.run(['f90wrap', '--version'], 
                         capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("ERROR: f90wrap not found. Please install it with:")
            print("  pip install f90wrap")
            sys.exit(1)

        # Run make to compile Fortran
        try:
            result = subprocess.run(['make', 'clean'], cwd=os.getcwd(), check=False)
            result = subprocess.run(['make', 'build'], cwd=os.getcwd(), check=True)
            
            # Copy compiled modules to src directory
            import shutil
            for ext in ['*.so', 'fieldline_tracer.py']:
                for f in Path('.').glob(ext):
                    shutil.copy(f, 'src/')
                    print(f"  Copied {f.name} to src/")
            
            print("\n" + "="*70)
            print("Fortran compilation completed successfully!")
            print("="*70 + "\n")
        except subprocess.CalledProcessError as e:
            print(f"\nERROR: Fortran compilation failed:")
            print(f"  {e}")
            sys.exit(1)

        # Don't call the default build_ext since we've handled it with Makefile
        # Just skip to parent run to complete any remaining setup
        super().run()


# Read long description from README
try:
    with open('README.md', 'r', encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "Magnetic field line tracer and VMEC interface for coil analysis"

# Read package metadata
package_name = 'coils2vmec'
try:
    with open(f'src/__init__.py', 'r') as f:
        for line in f:
            if line.startswith('__version__'):
                version = line.split('=')[1].strip().strip('"\'')
                break
        else:
            version = '0.2.0'
except FileNotFoundError:
    version = '0.2.0'

setup(
    name=package_name,
    version=version,
    description='Magnetic field line tracer and VMEC interface for coil analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Coils2VMEC Developers',
    author_email='your.email@example.com',
    url='https://github.com/yourusername/coils2vmec',
    package_dir={'': 'src'},
    packages=['coils2vmec'],
    cmdclass={'build_ext': FortranBuildExt},
    # Minimal ext_modules - actual compilation handled by Makefile
    ext_modules=[],
    zip_safe=False,
    include_package_data=True,
)
