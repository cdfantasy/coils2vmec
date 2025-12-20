# ==============================================================================
# Makefile for coils2vmec Python-Fortran interface
# ==============================================================================

# Compiler settings
FC = gfortran
FFLAGS = -O3 -fPIC -Wall -Wextra -fall-intrinsics
FFLAGS_DEBUG = -g -fPIC -Wall -Wextra -fcheck=all -fbacktrace
LDFLAGS = -shared

# Python settings
PYTHON = python3
F90WRAP = f90wrap
F2PY = f2py
PIP = $(PYTHON) -m pip

# Directories
FORTRAN_DIR = src/fortran
PYTHON_DIR = src/coils2vmec
BUILDDIR = build
TESTDIR = tests

# Absolute paths
FORTRAN_DIR_ABS := $(abspath $(FORTRAN_DIR))
PYTHON_DIR_ABS := $(abspath $(PYTHON_DIR))
BUILDDIR_ABS := $(abspath $(BUILDDIR))

# Source files
MAIN_MODULE = $(FORTRAN_DIR)/fieldline_tracer_module.f90
DEPENDENCIES = $(FORTRAN_DIR)/DLSODE.f
FORTRAN_SOURCES = $(MAIN_MODULE) $(DEPENDENCIES)
FORTRAN_OBJECTS = $(addsuffix .o, $(basename $(FORTRAN_SOURCES)))
FORTRAN_MODS = $(FORTRAN_DIR)/*.mod

# Python module name
MODULE_NAME = fieldline_tracer
SO_NAME = _$(MODULE_NAME)

# Generated files
WRAPPER_F90 = $(BUILDDIR)/f90wrap_fieldline_tracer_module.f90
PYTHON_WRAPPER = $(BUILDDIR)/$(MODULE_NAME).py
SO_FILE = $(BUILDDIR)/$(SO_NAME)*.so
TARGET_SO = $(PYTHON_DIR)/$(SO_NAME).so
TARGET_PY = $(PYTHON_DIR)/$(MODULE_NAME).py

# Build mode (release or debug)
BUILD_MODE ?= release
ifeq ($(BUILD_MODE),debug)
    FFLAGS = $(FFLAGS_DEBUG)
endif

# Version
VERSION = 0.1.0

# ==============================================================================
# Main targets
# ==============================================================================

.PHONY: all build fortran wrapper extension install clean clean_all help status check

all: build

build: fortran wrapper extension copy

# ==============================================================================
# Step 1: Compile Fortran source code
# ==============================================================================

fortran: $(FORTRAN_OBJECTS)
	@echo "✓ Step 1 complete: Fortran objects compiled"

$(FORTRAN_DIR)/%.o: $(FORTRAN_DIR)/%.f90
	@echo "Compiling Fortran 90: $<"
	$(FC) $(FFLAGS) -J$(FORTRAN_DIR) -c $< -o $@

$(FORTRAN_DIR)/%.o: $(FORTRAN_DIR)/%.f
	@echo "Compiling Fortran 77: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# Specific dependencies
$(FORTRAN_DIR)/fieldline_tracer_module.o: $(FORTRAN_DIR)/fieldline_tracer_module.f90
	@echo "Compiling main module: $<"
	$(FC) $(FFLAGS) -J$(FORTRAN_DIR) -c $< -o $@

$(FORTRAN_DIR)/DLSODE.o: $(FORTRAN_DIR)/DLSODE.f
	@echo "Compiling DLSODE library: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# ==============================================================================
# Step 2: Generate f90wrap wrapper code
# ==============================================================================

wrapper: $(WRAPPER_F90)

$(WRAPPER_F90): $(FORTRAN_SOURCES) kind_map.json
	@echo "✓ Step 2: Generating f90wrap interface..."
	@mkdir -p $(BUILDDIR)
	cd $(BUILDDIR) && $(F90WRAP) -m $(MODULE_NAME) \
		--kind-map ../kind_map.json \
		../$(MAIN_MODULE) \
		--verbose
	@echo "Fixing precision in generated wrappers..."
	@sed -i 's/real(4)/real(8)/g' $(BUILDDIR)/f90wrap_*.f90 2>/dev/null || \
	 sed -i '' 's/real(4)/real(8)/g' $(BUILDDIR)/f90wrap_*.f90
	@if [ -f "$(BUILDDIR)/$(MODULE_NAME).py" ]; then \
		echo "Fixing import statements in $(MODULE_NAME).py..."; \
		sed -i 's/^import _$(MODULE_NAME)$$/from . import _$(MODULE_NAME)/' $(BUILDDIR)/$(MODULE_NAME).py; \
		sed -i 's/^import _$(MODULE_NAME) as /from . import _$(MODULE_NAME) as /' $(BUILDDIR)/$(MODULE_NAME).py; \
		echo "✓ Import statements fixed"; \
	else \
		echo "Warning: $(MODULE_NAME).py not found for import fixing"; \
	fi
	@echo "✓ Wrapper files generated in $(BUILDDIR)"

# ==============================================================================
# Step 3: Build Python extension module (.so)
# ==============================================================================

extension: $(SO_FILE)

$(SO_FILE): $(FORTRAN_OBJECTS) $(WRAPPER_F90)
	@echo "✓ Step 3: Building Python extension module..."
	cd $(BUILDDIR) && $(F2PY) -c \
		--fcompiler=gnu95 \
		--f90exec=$(FC) \
		--f90flags="$(FFLAGS) -I$(FORTRAN_DIR_ABS)" \
		--opt="-O3" \
		-m $(SO_NAME) \
		f90wrap_fieldline_tracer_module.f90 \
		$(abspath $(FORTRAN_OBJECTS))
	@echo "✓ Extension module built in $(BUILDDIR)"

# ==============================================================================
# Step 4: Copy .so and .py to Python directory
# ==============================================================================

copy: $(TARGET_SO) $(TARGET_PY)

$(TARGET_SO): $(SO_FILE) | $(PYTHON_DIR)
	@echo "Copying .so file to Python directory..."
	cp $(BUILDDIR)/$(SO_NAME)*.so $(TARGET_SO)
	@echo "✓ Copied: $(TARGET_SO)"

$(TARGET_PY): $(PYTHON_WRAPPER) | $(PYTHON_DIR)
	@echo "Copying Python wrapper to Python directory..."
	cp $(PYTHON_WRAPPER) $(TARGET_PY)
	@echo "✓ Copied: $(TARGET_PY)"

# ==============================================================================
# Install via pip
# ==============================================================================

install:
	@echo "Installing package via pip..."
	$(PIP) install -e .
	@echo "✓ Package installed"

uninstall:
	@echo "Uninstalling package..."
	$(PIP) uninstall -y coils2vmec || true
	@echo "✓ Package uninstalled"

# ==============================================================================
# Testing
# ==============================================================================

test: build
	@echo "Testing module import..."
	$(PYTHON) -c "import sys; sys.path.insert(0, '$(PYTHON_DIR_ABS)'); import $(MODULE_NAME); print('✓ Successfully imported $(MODULE_NAME)')"

test_installed:
	@echo "Testing installed package..."
	$(PYTHON) -c "import coils2vmec.$(MODULE_NAME); print('✓ Successfully imported from installed package')"

# ==============================================================================
# Utility targets
# ==============================================================================

# Check dependencies
check:
	@echo "Checking dependencies..."
	@command -v $(FC) >/dev/null 2>&1 || { echo "✗ Error: $(FC) not found"; exit 1; }
	@command -v $(PYTHON) >/dev/null 2>&1 || { echo "✗ Error: $(PYTHON) not found"; exit 1; }
	@$(PYTHON) -c "import numpy" 2>/dev/null || { echo "✗ Error: numpy not installed"; exit 1; }
	@command -v $(F90WRAP) >/dev/null 2>&1 || { echo "✗ Error: f90wrap not found"; exit 1; }
	@command -v $(F2PY) >/dev/null 2>&1 || { echo "✗ Error: f2py not found"; exit 1; }
	@echo "✓ All dependencies found"

# Show build status
status:
	@echo "==================================="
	@echo "Build Status for coils2vmec"
	@echo "==================================="
	@echo ""
	@echo "1. Fortran objects:"
	@for obj in $(FORTRAN_OBJECTS); do \
		if [ -f $$obj ]; then \
			echo "  ✓ $$obj"; \
		else \
			echo "  ✗ $$obj (missing)"; \
		fi \
	done
	@echo ""
	@echo "2. Wrapper files:"
	@if [ -f $(WRAPPER_F90) ]; then echo "  ✓ $(WRAPPER_F90)"; else echo "  ✗ Wrapper .f90 not generated"; fi
	@if [ -f $(PYTHON_WRAPPER) ]; then echo "  ✓ $(PYTHON_WRAPPER)"; else echo "  ✗ Wrapper .py not generated"; fi
	@echo ""
	@echo "3. Extension module:"
	@if ls $(BUILDDIR)/$(SO_NAME)*.so >/dev/null 2>&1; then \
		echo "  ✓ Extension built in $(BUILDDIR)"; \
	else \
		echo "  ✗ Extension not built"; \
	fi
	@echo ""
	@echo "4. Python directory:"
	@if [ -f $(TARGET_SO) ]; then echo "  ✓ $(TARGET_SO)"; else echo "  ✗ .so not copied"; fi
	@if [ -f $(TARGET_PY) ]; then echo "  ✓ $(TARGET_PY)"; else echo "  ✗ .py not copied"; fi
	@echo ""

# Create necessary directories
$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

$(PYTHON_DIR):
	@mkdir -p $(PYTHON_DIR)

# Debug build
debug:
	$(MAKE) BUILD_MODE=debug build

# Rebuild from scratch
rebuild: clean build

# ==============================================================================
# Cleanup
# ==============================================================================

clean:
	@echo "Cleaning build files..."
	rm -f $(FORTRAN_DIR)/*.o $(FORTRAN_DIR)/*.mod
	rm -rf $(BUILDDIR)
	rm -f $(PYTHON_DIR)/*.so
	rm -f $(PYTHON_DIR)/$(MODULE_NAME).py
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	@echo "✓ Build files cleaned"

clean_all: clean
	@echo "Cleaning all generated files..."
	rm -rf dist *.egg-info
	rm -rf .pytest_cache
	@echo "✓ All files cleaned"

distclean: clean_all uninstall
	@echo "✓ Full clean complete"

# ==============================================================================
# Help
# ==============================================================================

help:
	@echo "======================================"
	@echo "coils2vmec Makefile"
	@echo "======================================"
	@echo ""
	@echo "Project structure:"
	@echo "  Fortran source:  $(FORTRAN_DIR)/"
	@echo "  Python package:  $(PYTHON_DIR)/"
	@echo "  Build directory: $(BUILDDIR)/"
	@echo ""
	@echo "Build workflow:"
	@echo "  1. make          - Build extension (Fortran + wrapper + .so)"
	@echo "  2. pip install . - Register to system (run separately)"
	@echo ""
	@echo "  Details:"
	@echo "    Step 1: Compile Fortran → .o files"
	@echo "    Step 2: Generate f90wrap wrapper → .f90, .py"
	@echo "    Step 3: Build extension → .so file"
	@echo "    Step 4: Copy to Python directory"
	@echo ""
	@echo "Available targets:"
	@echo "  make               - Build extension (default)"
	@echo "  make build         - Complete build process"
	@echo "  make fortran       - Compile Fortran code only"
	@echo "  make wrapper       - Generate f90wrap wrappers"
	@echo "  make extension     - Build .so extension"
	@echo "  make copy          - Copy files to Python dir"
	@echo "  make install       - Install via pip (run after make)"
	@echo "  make uninstall     - Uninstall package"
	@echo ""
	@echo "Testing:"
	@echo "  make test          - Test local build"
	@echo "  make test_installed - Test installed package"
	@echo ""
	@echo "Utilities:"
	@echo "  make check         - Check dependencies"
	@echo "  make status        - Show build status"
	@echo "  make debug         - Debug build"
	@echo "  make rebuild       - Clean and rebuild"
	@echo ""
	@echo "Cleanup:"
	@echo "  make clean         - Clean build files"
	@echo "  make clean_all     - Clean all generated files"
	@echo "  make distclean     - Full clean + uninstall"
	@echo ""
	@echo "Options:"
	@echo "  BUILD_MODE=debug   - Build with debug flags"
	@echo ""
	@echo "Compiler: $(FC) $(FFLAGS)"
	@echo "Python:   $(PYTHON)"

.DEFAULT_GOAL := all
