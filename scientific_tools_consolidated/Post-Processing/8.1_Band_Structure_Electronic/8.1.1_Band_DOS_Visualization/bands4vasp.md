# bands4vasp

## Official Resources
- **Homepage**: https://github.com/QuantumMaterialsModelling/bands4vasp
- **GitHub**: https://github.com/QuantumMaterialsModelling/bands4vasp
- **Publication**: J. Phys. Chem. C 2021, 125, 21, 11714-11724
- **License**: GPL v3

## Overview
bands4vasp is a post-processing package for the analysis of unfolded eigenstates in VASP. It provides comprehensive tools for band structures, 2D and 3D Fermi surfaces, Fermi vectors, and spectral functions from VASP calculations with the unfolding patch.

**Scientific domain**: Band unfolding, Fermi surface analysis, spectral functions
**Target user community**: VASP users studying supercells, alloys, and disordered systems

## Theoretical Background
bands4vasp processes unfolded VASP data based on:
- Spectral weight from supercell unfolding
- Fermi surface extraction at E = E_F
- Fermi vector determination in k-space
- Energy-resolved spectral function A(k,E)

## Capabilities (CRITICAL)
- **Unfolded Bands**: Supercell band structure unfolding visualization
- **Fermi Surfaces**: 2D and 3D Fermi surface visualization
- **Spectral Functions**: Energy-resolved spectral analysis A(k,E)
- **Fermi Vectors**: k-space Fermi vector extraction
- **Spin Textures**: Spin-resolved analysis
- **Publication Quality**: High-quality matplotlib output

## Key Strengths

### VASP Integration:
- Native VASP 6.2.1+ support
- Unfolding patch integration
- Direct output file reading

### Visualization:
- 2D band structure plots
- 3D Fermi surface rendering
- Spectral function heatmaps
- Publication-ready figures

### Analysis:
- Fermi vector extraction
- Spectral weight analysis
- Spin texture mapping

## Requirements
- VASP 6.2.1+ with unfolding patch
- Python 3.x
- NumPy, Matplotlib

## Inputs & Outputs
- **Input formats**:
  - VASP unfolding output files
  - EIGENVAL with unfolding data
  - PROCAR for projections
  
- **Output data types**:
  - Band structure plots
  - Fermi surface visualizations
  - Spectral function images
  - Numerical data exports

## Installation
```bash
git clone https://github.com/QuantumMaterialsModelling/bands4vasp.git
cd bands4vasp
pip install -e .
```

## Usage Examples
```python
from bands4vasp import BandStructure, FermiSurface

# Plot unfolded band structure
bs = BandStructure.from_vasp("path/to/vasp/calc")
bs.plot()

# Plot Fermi surface
fs = FermiSurface.from_vasp("path/to/vasp/calc")
fs.plot_3d()
```

## Performance Characteristics
- **Speed**: Efficient post-processing
- **Memory**: Handles large unfolding datasets
- **Visualization**: High-quality matplotlib output

## Limitations & Known Constraints
- **VASP-specific**: Requires VASP with unfolding patch
- **Version requirement**: VASP 6.2.1 or later
- **Patch needed**: Unfolding functionality requires patched VASP

## Comparison with Other Tools
- **vs easyunfold**: bands4vasp uses VASP native unfolding
- **vs BandUP**: Different unfolding approach
- **Unique strength**: Native VASP unfolding integration, Fermi surfaces

## Application Areas
- Alloy band structures
- Defect electronic structure
- Disordered system analysis
- Fermi surface studies
- Spectral function analysis

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/QuantumMaterialsModelling/bands4vasp
2. J. Phys. Chem. C 2021, 125, 21, 11714-11724

**Confidence**: VERIFIED - Published in peer-reviewed journal

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Developer: Quantum Materials Modelling group
- Academic citations: Published in J. Phys. Chem. C
