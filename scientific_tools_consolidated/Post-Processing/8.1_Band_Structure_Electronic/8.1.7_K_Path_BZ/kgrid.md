# kgrid

## Overview
**kgrid** is a Python tool for calculating the required k-point density from input geometry for periodic quantum chemistry calculations.

## Official Resources
- **GitHub**: https://github.com/WMD-group/kgrid
- **PyPI**: https://pypi.org/project/kgrid/

## Capabilities
- **K-point Generation**: Calculate optimal k-grids
- **Density-based**: Specify k-point density
- **Multi-code Support**: VASP, CASTEP output formats

## Key Features
- Length cutoff-based k-point selection
- KSPACING output for VASP
- MP SPACING for CASTEP
- Command-line interface

## Usage
```bash
kgrid POSCAR 25  # 25 Å cutoff
kgrid --castep POSCAR 25  # CASTEP format
```

## Installation
```bash
pip install kgrid
```

## Verification & Sources
- **Status**: ✅ VERIFIED
- **Confidence**: VERIFIED
- **Developer**: WMD Group
