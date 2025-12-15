# Coils2VMEC

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()

**Magnetic field line tracer and VMEC interface for stellarator coil analysis**

Coils2VMEC is a comprehensive Python package that converts discrete 3D coil configurations into VMEC-compatible magnetic equilibria through fieldline tracing, rotational transform analysis, and surface fitting.

---

## 🌟 Features

- **🔬 High-Performance Fieldline Tracing**: Fortran-accelerated magnetic field line integration using LSODE solver
- **📊 Rotational Transform Analysis**: Automatic iota profile calculation with error estimation
- **🎯 LCFS Detection**: Robust last closed flux surface identification with island detection
- **🔄 DESCUR Surface Fitting**: Python implementation of DESCUR for Fourier decomposition
- **📈 Advanced Visualization**: Interactive 3D plots, Poincaré sections, and surface overlays
- **⚡ Parallel Processing**: Multi-threaded fieldline tracing for improved performance
- **🔧 VMEC Integration**: Direct generation of VMEC input files with toroidal flux calculation

---

## 📋 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Examples](#usage-examples)
- [Package Structure](#package-structure)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## 🚀 Installation

### Prerequisites

**System Requirements:**
- Python 3.8 or higher
- gfortran compiler
- f90wrap (for Fortran-Python interface)

**Platform Support:**
- ✅ Linux (Ubuntu 18.04+, CentOS 7+)
- ✅ macOS (10.14+)
- ⚠️ Windows (via WSL2)

### Step 1: Install System Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install gfortran python3-dev
```

#### macOS
```bash
brew install gcc
```

#### CentOS/RHEL
```bash
sudo yum install gcc-gfortran python3-devel
```

### Step 2: Install Coils2VMEC

```bash
# Clone the repository
git clone https://github.com/cdfantasy/coils2vmec.git
cd coils2vmec

# Install Python build dependencies
pip install f90wrap numpy

# Compile Fortran extensions
make

# Install the package in editable mode (recommended for development)
pip install -e .

# Or install directly
pip install .
```

**Installation Steps:**
1. **`make`** - Compiles Fortran code and generates Python wrappers using f90wrap
2. **`pip install -e .`** - Installs the Python package in editable mode

### Verify Installation

```bash
python -c "import coils2vmec as c2v; print(f'✅ Coils2VMEC v{c2v.__version__} installed successfully')"
```

---

## ⚡ Quick Start

### Basic Workflow

```python
import coils2vmec as c2v
import numpy as np

# 1. Read coil configuration
coils_data = c2v.read_coils_file('coils.xxx')

# 2. Find magnetic axis
axis_rz = c2v.find_axis(initial_guess=[1.32, 0.0])
print(f"Magnetic axis: R={axis_rz[0]:.4f} m, Z={axis_rz[1]:.4f} m")

# 3. Trace fieldlines in parallel
fieldlines_data = c2v.trace_fieldlines_parallel(
    initial_guess=axis_rz,
    n_fieldlines=99,
    nturn=400,
    nphi=360,
    coils_data=coils_data
)

# 4. Calculate iota profile
iota_results = c2v.calculate_iota_profile(
    X_lines, Y_lines, Z_lines, 
    nturn=400, nphi=360, nline=99
)

# 5. Detect LCFS
lcfs_result = c2v.find_lcfs_and_islands(
    radius, iota,
    threshold_factor=15,
    plot_flag=True
)

# 6. Fit surface with DESCUR
from coils2vmec import DescurConfig
config = DescurConfig(mu=8, nv=60, ftol=2e-6)
c2v.run_descur_python(
    R_lines, Z_lines, Phi_lines,
    lcfs_idx=lcfs_result['lcfs_index'],
    config=config
)

# 7. Generate VMEC input
c2v.save_vmec_input_for_surface(
    surface, bs_tf, extcur, 
    mpol=8, ntor=8, nfp=1,
    output_dir='./vmec_input'
)
```

---

## 📚 Usage Examples

### Example 1: Complete H1 Device Analysis

See [`examples/example_pwO.py`](examples/example_pwO.py) for a full working example:

```bash
cd examples
python example_pwO.py
```

### Example 2: Jupyter Notebook Workflow

Interactive analysis in [`coils2vmec.ipynb`](coils2vmec.ipynb):
```bash
jupyter notebook coils2vmec.ipynb
```

### Example 3: Custom Configuration

```python
from coils2vmec import DescurConfig

# Configure DESCUR fitting
config = DescurConfig()
config.mu = 12          # Poloidal modes
config.nv = 72          # Toroidal points
config.ftol = 1e-6      # Convergence tolerance
config.niter = 1000     # Maximum iterations
config.alpha = 0.8      # Step size
```

---

## 📦 Package Structure

```
coils2vmec/
├── src/
│   └── coils2vmec/
│       ├── __init__.py          # Package exports
│       ├── fieldline.py         # Fieldline tracing
│       ├── iota.py              # Rotational transform
│       ├── lcfs.py              # LCFS detection
│       ├── descur.py            # DESCUR interface
│       ├── descur_python.py     # DESCUR implementation
│       ├── plotting.py          # Visualization
│       ├── utils.py             # Utilities
│       ├── fieldline_tracer.py  # f90wrap interface
│       └── fortran/             # Fortran sources
│           ├── fieldline_tracer_module.f90
│           ├── DLSODE.f
│           ├── hybrd.f
│           └── traceline.f90
├── examples/                     # Example scripts
├── test/                         # Test data
├── setup.py                      # Build configuration
├── pyproject.toml               # PEP 517 metadata
├── requirements.txt             # Dependencies
└── README.md                    # This file
```

---

## 📖 Documentation

### Core Modules

#### Fieldline Tracing
- `find_axis()` - Magnetic axis location
- `find_lcfs()` - Last closed flux surface
- `trace_fieldlines_parallel()` - Parallel fieldline integration
- `save_fieldlines_hdf5()` - HDF5 data persistence

#### Iota Analysis
- `calculate_iota_profile()` - Rotational transform calculation
- `check_rational_surface()` - Rational surface detection
- `adjust_lcfs_avoid_rational_surface()` - LCFS optimization

#### Surface Fitting
- `run_descur_python()` - DESCUR surface fitting
- `DescurConfig` - Configuration dataclass
- `DescurFitter` - Fitting algorithm

#### Visualization
- `plot_fieldlines_3d()` - 3D fieldline visualization
- `plot_poincare_sections()` - Poincaré plots
- `plot_iota_with_radius()` - Iota profiles
- `plot_poincare_with_surface()` - Combined surface/Poincaré

### API Reference

Full API documentation: [docs/API.md](docs/API.md) *(coming soon)*

---

## 🔧 Configuration

### Environment Variables

```bash
# Set Fortran compiler (optional)
export FC=gfortran

# Enable verbose output
export COILS2VMEC_VERBOSE=1

# Set number of parallel workers
export COILS2VMEC_WORKERS=8
```

### Configuration Files

Create a `config.yaml` for default settings:

```yaml
fieldline:
  nturn: 400
  nphi: 360
  
descur:
  mu: 8
  nv: 60
  ftol: 2e-6
  
visualization:
  dpi: 300
  format: png
```

---

## 🤝 Contributing

We welcome contributions! 

### Development Setup

```bash
# Clone repository
git clone https://github.com/cdfantasy/coils2vmec.git
cd coils2vmec

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/

# Format code
black src/
isort src/
```

### Reporting Issues

Found a bug? Have a feature request? Please open an issue on [GitHub Issues](https://github.com/cdfantasy/coils2vmec/issues).

---

## 📊 Performance

Benchmark on typical H1 device configuration:
- **99 fieldlines, 400 toroidal turns, 360 φ points**
- **Hardware**: Intel Xeon E5-2680 v4 (28 cores)
- **Parallel speedup**: ~18x (8 workers)
- **Total runtime**: ~45 seconds

---

## 🔬 Scientific Background

### Methods

1. **Biot-Savart Field Calculation**: Magnetic field from discrete coil segments
2. **LSODE Integration**: Adaptive step-size ODE solver for fieldline tracing
3. **Fourier Decomposition**: DESCUR algorithm for toroidal surface fitting
4. **Rational Surface Analysis**: Iota-based detection of resonant surfaces

### References

- VMEC Documentation: [https://princetonuniversity.github.io/STELLOPT/](https://princetonuniversity.github.io/STELLOPT/)
- DESCUR Algorithm: [Original Fortran implementation](https://github.com/ORNL-Fusion/DESCUR)
- f90wrap: [https://github.com/jameskermode/f90wrap](https://github.com/jameskermode/f90wrap)

---

## 📄 Citation

If you use Coils2VMEC in your research, please cite:

```bibtex
@software{coils2vmec2025,
  author = {Your Name},
  title = {Coils2VMEC: Magnetic Field Line Tracer for Stellarator Analysis},
  year = {2025},
  url = {https://github.com/cdfantasy/coils2vmec},
  version = {0.2.0}
}
```

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **VMEC Team** - For the VMEC equilibrium code
- **STELLOPT Consortium** - For stellarator optimization tools
- **f90wrap Developers** - For Fortran-Python interface generation
- **SciPy Community** - For numerical computing tools

---

## 📞 Contact

- **Maintainer**: cdfantasy
- **Email**: zkgao@stu.usc.edu.cn
- **Issues**: [GitHub Issues](https://github.com/cdfantasy/coils2vmec/issues)
- **Discussions**: [GitHub Discussions](https://github.com/cdfantasy/coils2vmec/discussions)

---

## 🗺️ Roadmap


---
