# postqe

## Official Resources
- **Homepage**: https://github.com/QEF/postqe
- **GitHub**: https://github.com/QEF/postqe
- **Documentation**: https://postqe.readthedocs.io/
- **PyPI**: https://pypi.org/project/postqe/
- **License**: GPL v2

## Overview
postqe is a Python package for post-processing calculations from Quantum ESPRESSO, developed by the Quantum ESPRESSO Foundation. It provides functions to analyze and visualize results including charge density, potentials, electronic structure, and equation of state fitting.

**Scientific domain**: Quantum ESPRESSO post-processing, electronic structure visualization
**Target user community**: Quantum ESPRESSO users requiring Python-based analysis

## Theoretical Background
postqe processes Quantum ESPRESSO output to analyze:
- Charge density ρ(r) on 1D/2D sections
- Electrostatic potentials (bare, Hartree, total)
- Kohn-Sham eigenvalues and band structures
- Equation of state E(V) fitting

## Capabilities (CRITICAL)
- **Charge Density**: Plot on 1D/2D sections through crystal
- **Potentials**: Bare/Hartree/total potential visualization
- **EOS Fitting**: Murnaghan, Birch-Murnaghan equations of state
- **Band Structure**: Electronic band plotting
- **DOS**: Density of states analysis
- **XML Parsing**: Read QE XML output files

## Key Strengths

### Native QE Support:
- Direct XML file parsing
- All QE output types supported
- Official QEF development

### Visualization:
- 1D line plots through structures
- 2D contour plots
- Matplotlib integration
- Publication-quality figures

### EOS Analysis:
- Multiple EOS models
- Bulk modulus extraction
- Equilibrium volume fitting

## Inputs & Outputs
- **Input formats**:
  - QE XML output files
  - Charge density files
  - Potential files
  - Band structure data
  
- **Output data types**:
  - Matplotlib figures
  - Numerical data arrays
  - Fitted parameters

## Installation
```bash
pip install postqe
```

## Usage Examples
```python
from postqe import get_charge, plot_charge_1d

# Read charge density
charge = get_charge("charge-density.dat")

# Plot 1D section
plot_charge_1d(charge, [0,0,0], [1,0,0])

# EOS fitting
from postqe import fit_eos
volumes, energies = [...], [...]
params = fit_eos(volumes, energies, eos='birch-murnaghan')
```

## Performance Characteristics
- **Speed**: Efficient file parsing
- **Memory**: Handles large charge density files
- **Visualization**: Fast matplotlib rendering

## Limitations & Known Constraints
- **QE-specific**: Only for Quantum ESPRESSO output
- **XML format**: Requires QE XML output enabled
- **Limited features**: Focused on basic post-processing

## Comparison with Other Tools
- **vs abipy**: postqe for QE, abipy for ABINIT
- **vs pymatgen**: postqe more QE-specific
- **Unique strength**: Official QEF tool, EOS fitting

## Application Areas
- Charge density analysis
- Potential visualization
- Equation of state studies
- Band structure plotting
- Basic QE post-processing

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/QEF/postqe

**Confidence**: VERIFIED - Official QEF tool

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL v2)
- Developer: QEF (Quantum ESPRESSO Foundation)
