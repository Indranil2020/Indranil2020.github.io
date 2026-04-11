# GAP / QUIP

## Official Resources
- Source Repository: https://github.com/libAtoms/QUIP
- Documentation: https://libatoms.github.io/QUIP/
- License: Open source (GPL-2.0)

## Overview
**GAP** (Gaussian Approximation Potential) with **QUIP** (Quantum Interatomic Potential) is the original MLIP framework using SOAP descriptors and sparse Gaussian processes. It pioneered the field of ML interatomic potentials and remains widely used for element-specific high-accuracy potentials.

**Scientific domain**: Gaussian process MLIP with SOAP descriptors  
**Target user community**: Researchers fitting high-accuracy GP potentials for specific elements/systems

## Theoretical Methods
- Gaussian Approximation Potential (GAP)
- SOAP (Smooth Overlap of Atomic Positions) descriptors
- Sparse Gaussian process regression
- 2-body, 3-body, and SOAP terms
- Fortran core with Python interface (quippy)

## Capabilities (CRITICAL)
- SOAP-GAP fitting
- Multi-body terms (2b, 3b, SOAP)
- LAMMPS integration via quippy
- ASE interface
- Well-tested on many systems

**Sources**: GitHub repository, Int. J. Quantum Chem. 115, 1051 (2015)

## Key Strengths

### Pioneering:
- First widely-used MLIP framework
- SOAP descriptor standard
- Well-validated on many systems
- Published potentials library

### Accuracy:
- Near-DFT accuracy
- Well-calibrated uncertainty
- Systematic improvement
- Physical constraints

### Integration:
- LAMMPS via quippy
- ASE calculator
- Fortran performance
- Python flexibility

## Inputs & Outputs
- **Input formats**: Extended XYZ, VASP, CASTEP
- **Output data types**: GAP potentials, energies, forces, virials

## Interfaces & Ecosystem
- **LAMMPS**: MD engine
- **ASE**: Calculator (quippy)
- **Fortran**: Core
- **Python**: Interface

## Performance Characteristics
- **Speed**: Moderate (GP evaluation)
- **Accuracy**: Near-DFT
- **System size**: 100-10000 atoms
- **Automation**: Semi-automated

## Computational Cost
- **Fitting**: Hours
- **MD**: 10-100x faster than DFT

## Limitations & Known Constraints
- **GPL license**: Copyleft
- **Fortran build**: Complex compilation
- **GP scaling**: O(N²) training
- **Element-specific**: Not universal out-of-box

## Comparison with Other Codes
- **vs MACE**: GAP is GP, MACE is NN
- **vs FLARE**: GAP is batch, FLARE is on-the-fly
- **vs SNAP**: GAP is GP, SNAP is linear
- **Unique strength**: Pioneering MLIP framework with SOAP descriptors and sparse GP, well-validated across many systems

## Application Areas

### High-Accuracy Potentials:
- Carbon (GAP-20)
- Silicon, germanium
- Tungsten, iron
- Water, aqueous systems

### Materials MD:
- Phase transitions
- Amorphous materials
- Defect properties
- Thermal transport

## Best Practices
- Use SOAP + 2b + 3b terms
- Start from published potentials
- Validate with phonons and elastic constants
- Use quippy for LAMMPS integration

## Community and Support
- Open source (GPL-2.0)
- libAtoms community
- Cambridge/Oxford maintained
- Comprehensive documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/libAtoms/QUIP

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Pioneering MLIP framework with SOAP descriptors and sparse GP
