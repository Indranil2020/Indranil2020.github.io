# banduppy

## Official Resources
- **Homepage**: https://github.com/band-unfolding/banduppy
- **GitHub**: https://github.com/band-unfolding/banduppy
- **Original BandUP**: https://github.com/band-unfolding/bandup
- **PyPI**: https://pypi.org/project/banduppy/
- **License**: GPL v3

## Overview
banduppy is a Python implementation of the BandUP code for band structure unfolding from supercell calculations. It provides modern support for Quantum ESPRESSO and other DFT codes, using the irrep library for wavefunction reading. banduppy enables extraction of effective band structures from supercell calculations containing defects, alloys, or interfaces.

**Scientific domain**: Band structure unfolding, supercell calculations, defect physics
**Target user community**: Researchers studying defects, alloys, and disordered systems using supercell methods

## Theoretical Background
banduppy implements band unfolding based on:
- Spectral weight: W(k,E) = Σ_n |⟨ψ_nK|e^{ik·r}⟩|² δ(E - E_nK)
- Maps supercell K-points to primitive cell k-points
- Uses irreducible representations for symmetry analysis
- Based on Popescu & Zunger methodology

## Capabilities (CRITICAL)
- **Band Unfolding**: Supercell to primitive cell unfolding
- **Spectral Weights**: Calculate and visualize unfolding weights
- **Multi-code Support**: QE, VASP, ABINIT, CASTEP
- **Wavefunction Reading**: Uses irrep library routines
- **Spectral Functions**: Energy-resolved spectral analysis
- **Effective Bands**: Extract effective band structure

## Key Strengths

### Python Implementation:
- Pure Python, easy installation
- Modern Python 3 support
- Integration with scientific Python stack
- Scriptable workflows

### Multi-Code Support:
- Quantum ESPRESSO (primary)
- VASP
- ABINIT
- CASTEP
- Extensible architecture

### irrep Integration:
- Symmetry-aware analysis
- Irreducible representations
- Robust wavefunction reading

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO output
  - VASP WAVECAR
  - ABINIT WFK files
  - CASTEP output
  
- **Output data types**:
  - Spectral weights
  - Unfolded band structures
  - Spectral functions

## Installation
```bash
pip install banduppy
```

## Usage Examples
```python
from banduppy import UnfoldingPath, Unfolding

# Define unfolding path
path = UnfoldingPath(
    supercell_matrix=[[2,0,0],[0,2,0],[0,0,2]],
    kpath_primitive=[[[0,0,0],[0.5,0,0]]]
)

# Perform unfolding
unfold = Unfolding(path, wavefunction_file="wfc.dat")
spectral_weights = unfold.get_spectral_weights()
```

## Performance Characteristics
- **Speed**: Efficient Python implementation
- **Memory**: Handles large supercells
- **Accuracy**: Proper symmetry handling via irrep

## Limitations & Known Constraints
- **Wavefunction required**: Needs wavefunction output
- **Memory**: Large supercells need significant RAM
- **Code-specific**: Different interfaces per code

## Comparison with Other Tools
- **vs BandUP (Fortran)**: banduppy is Python, easier installation
- **vs easyunfold**: Different approaches, both effective
- **Unique strength**: irrep integration, multi-code support

## Application Areas
- Defect band structures
- Alloy electronic structure
- Interface/heterostructure bands
- Disordered system analysis

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/band-unfolding/banduppy
2. Original BandUP: https://github.com/band-unfolding/bandup

**Confidence**: VERIFIED - Python version of established BandUP code

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Active development: Maintained
