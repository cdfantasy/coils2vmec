#!/usr/bin/env python3
"""
Setup script for coils2vmec package

Installation:
    1. Compile Fortran extensions: make
    2. Install Python package: pip install -e .
"""

from pathlib import Path
from setuptools import setup

# Read long description from README
try:
    with open('README.md', 'r', encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "Magnetic field line tracer and VMEC interface for coil analysis"

# Read package metadata
package_name = 'coils2vmec'
try:
    with open('src/coils2vmec/__init__.py', 'r') as f:
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
    author_email='zkgao@stu.usc.edu.cn',
    url='https://github.com/cdfantasy/coils2vmec',
    package_dir={'': 'src'},
    packages=['coils2vmec'],
    zip_safe=False,
    include_package_data=True,
)
