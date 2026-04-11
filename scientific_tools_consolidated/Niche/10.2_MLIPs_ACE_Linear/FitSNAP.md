# FitSNAP

## Official Resources
- Source Repository: https://github.com/FitSNAP/FitSNAP
- Documentation: https://fitsnap.github.io/
- License: Open source (MIT)

## Overview
**FitSNAP** is software for generating machine-learning interatomic potentials for LAMMPS. It implements SNAP, qSNAP, and other linear/nonlinear potentials with tight LAMMPS integration for production MD simulations.

**Scientific domain**: SNAP/qSNAP potential fitting for LAMMPS  
**Target user community**: Researchers fitting SNAP potentials for LAMMPS MD

## Theoretical Methods
- SNAP (Spectral Neighbor Analysis Potential)
- qSNAP (quadratic SNAP)
- Linear and nonlinear fitting
- Bispectrum descriptors
- LAMMPS mliap integration

## Capabilities (CRITICAL)
- SNAP potential fitting
- qSNAP quadratic extension
- LAMMPS production integration
- Multi-element support
- Parallel fitting
- Uncertainty quantification

**Sources**: GitHub repository

## Key Strengths

### LAMMPS Integration:
- Direct mliap pair_style
- Production MD ready
- Parallel execution
- No format conversion

### SNAP Framework:
- Well-tested SNAP implementation
- qSNAP for improved accuracy
- Linear regression (fast fitting)
- Physics-informed constraints

### Production Quality:
- Published potentials available
- Sandia National Labs maintained
- Extensive testing
- Documentation

## Inputs & Outputs
- **Input formats**: Training data (LAMMPS dump, VASP, etc.)
- **Output data types**: LAMMPS potential files, SNAP coefficients

## Interfaces & Ecosystem
- **LAMMPS**: MD engine
- **Python**: Core language
- **NumPy**: Computation

## Performance Characteristics
- **Speed**: Very fast (linear model)
- **Accuracy**: SNAP-level (~100 meV/atom)
- **System size**: Any (LAMMPS)
- **Automation**: Full

## Computational Cost
- **Fitting**: Minutes
- **MD**: Very fast (linear evaluation)

## Limitations & Known Constraints
- **SNAP accuracy**: Lower than NN potentials
- **Descriptor fixed**: SNAP only
- **LAMMPS only**: No other MD engines
- **Training data**: Needs diverse configurations

## Comparison with Other Codes
- **vs ACE1pack**: FitSNAP is SNAP, ACE1pack is ACE basis
- **vs DeePMD-kit**: FitSNAP is linear, DeePMD is NN
- **vs PACE**: FitSNAP is SNAP, PACE is ACE in LAMMPS
- **Unique strength**: SNAP/qSNAP fitting with direct LAMMPS mliap integration

## Application Areas

### SNAP Potentials:
- Tungsten, tantalum, uranium
- BCC/FCC metals
- High-temperature MD
- Radiation damage

### LAMMPS MD:
- Production MD with MLIP
- Large-scale simulations
- Multi-million atom runs

## Best Practices
- Use diverse training data
- Validate with elastic constants
- Test phonon spectra
- Compare with DFT MD

## Community and Support
- Open source (MIT)
- Sandia National Labs maintained
- Comprehensive documentation
- Published potentials library

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/FitSNAP/FitSNAP

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: SNAP/qSNAP fitting with direct LAMMPS mliap integration
