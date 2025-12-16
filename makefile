# ==============================================================================
# Makefile for fieldline_tracer Python interface using f90wrap
# ==============================================================================

# Compiler settings
FC = gfortran
FFLAGS = -O2 -fPIC -Wall -Wextra -fall-intrinsics
LDFLAGS = -shared

# Python settings
PYTHON = python3
F90WRAP = f90wrap
F2PY = f2py

# Directories
SRCDIR = src/coils2vmec/fortran
BUILDDIR = build
WRAPDIR = f90wrap_generated
TESTDIR = tests
DISTDIR = dist

# Source files
MAIN_MODULE = $(SRCDIR)/fieldline_tracer_module.f90
DEPENDENCIES = $(SRCDIR)/DLSODE.f
FORTRAN_SOURCES = $(MAIN_MODULE) $(DEPENDENCIES)
FORTRAN_OBJECTS = $(addsuffix .o, $(basename $(FORTRAN_SOURCES)))

# Python module name
MODULE_NAME = fieldline_tracer
PYTHON_MODULE = _$(MODULE_NAME)*.so

# Generated files
WRAPPER_F90 = f90wrap_fieldline_tracer_module.f90
PYTHON_WRAPPER = $(MODULE_NAME).py

# Version info
VERSION = 0.1.0

# Detect OS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    SHLIB_EXT = .so
endif
ifeq ($(UNAME_S),Darwin)
    SHLIB_EXT = .dylib
endif

# ==============================================================================
# Main targets
# ==============================================================================

.PHONY: all build install test clean clean_all help

all: build

build: $(PYTHON_MODULE)

# ==============================================================================
# Python module build
# ==============================================================================

# Generate f90wrap interface
$(WRAPPER_F90): $(FORTRAN_SOURCES) kind_map.json
	@echo "Step 1: Generating f90wrap interface..."
	$(F90WRAP) -m $(MODULE_NAME) \
		--kind-map kind_map.json \
		$(MAIN_MODULE) \
		--verbose
	@echo "Fixing precision in generated wrappers..."
	@sed -i 's/real(4)/real(8)/g' f90wrap_*.f90 || sed -i '' 's/real(4)/real(8)/g' f90wrap_*.f90
	@echo "f90wrap files generated"

# Compile Python extension module
$(PYTHON_MODULE): $(FORTRAN_OBJECTS) $(WRAPPER_F90)
	@echo "Step 2: Building Python extension module..."
	$(F2PY) -c \
		--fcompiler=gnu95 \
		--f90exec=$(FC) \
		--f90flags="$(FFLAGS) -I. -I$(SRCDIR)" \
		--opt="-O3" \
		-m _$(MODULE_NAME) \
		$(WRAPPER_F90) \
		$(FORTRAN_OBJECTS) \
		--build-dir $(BUILDDIR)
	@echo "Step 3: Python module built successfully!"

# ==============================================================================
# Fortran compilation
# ==============================================================================

# Pattern rules
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# Dependencies
$(SRCDIR)/fieldline_tracer_module.o: $(SRCDIR)/fieldline_tracer_module.f90
	$(FC) $(FFLAGS) -J$(SRCDIR) -c $< -o $@

$(SRCDIR)/DLSODE.o: $(SRCDIR)/DLSODE.f
	$(FC) $(FFLAGS) -c $< -o $@

# Compile Fortran objects separately
fortran_objects: $(FORTRAN_OBJECTS)
	@echo "Fortran objects compiled"

# ==============================================================================
# Testing
# ==============================================================================

test: build
	@echo "Testing module import..."
	$(PYTHON) -c "import $(MODULE_NAME); print('✓ Successfully imported $(MODULE_NAME)')"
	@echo "Testing with sample data..."
	$(PYTHON) test_import.py

test_fortran: fortran_objects
	@echo "Testing Fortran compilation..."
	$(FC) $(FFLAGS) -o test_fortran test_fortran.f90 $(FORTRAN_OBJECTS)
	@echo "Fortran test program compiled"

# ==============================================================================
# Installation
# ==============================================================================

install: setup.py
	$(PYTHON) -m pip install .

develop:
	$(PYTHON) -m pip install -e .

# ==============================================================================
# Distribution
# ==============================================================================

sdist: setup.py
	$(PYTHON) setup.py sdist

bdist_wheel: setup.py
	$(PYTHON) setup.py bdist_wheel

# ==============================================================================
# Cleanup
# ==============================================================================

clean:
	@echo "Cleaning build files..."
	rm -f *.o *.mod *.so *.dylib
	rm -f f90wrap_*.f90 $(MODULE_NAME).py .f2py_f2cmap
	rm -rf $(BUILDDIR) $(WRAPDIR) __pycache__
	rm -rf coils2vmec/  # Remove f2py generated package directory
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

clean_all: clean
	@echo "Cleaning all generated files..."
	rm -f $(PYTHON_MODULE)
	rm -rf $(DISTDIR) *.egg-info
	rm -f setup.py

distclean: clean_all
	@echo "Full clean complete"

# ==============================================================================
# Help
# ==============================================================================

help:
	@echo "Field Line Tracer - Python Interface Makefile"
	@echo "============================================="
	@echo ""
	@echo "Available commands:"
	@echo "  make               - Build Python module (default)"
	@echo "  make build         - Build Python module"
	@echo "  make fortran_objects - Compile Fortran objects only"
	@echo "  make test          - Test the built module"
	@echo "  make install       - Install via pip"
	@echo "  make develop       - Install in development mode"
	@echo "  make clean         - Clean build files"
	@echo "  make clean_all     - Clean all generated files"
	@echo "  make help          - Show this help"
	@echo ""
	@echo "Source files:"
	@echo "  Fortran: $(FORTRAN_SOURCES)"
	@echo "  Python module: $(MODULE_NAME)"
	@echo ""
	@echo "Compiler: $(FC) $(FFLAGS)"

.DEFAULT_GOAL := all