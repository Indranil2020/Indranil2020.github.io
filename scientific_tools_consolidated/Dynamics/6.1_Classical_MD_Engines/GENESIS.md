# GENESIS

## Official Resources
- Homepage: https://www.r-ccs.riken.jp/labs/cbrt/
- Documentation: https://www.r-ccs.riken.jp/labs/cbrt/genesis/
- Source Repository: https://github.com/genesis-release-r-ccs/genesis
- License: LGPL

## Overview
GENESIS (GENeralized-Ensemble SImulation System) is a high-performance molecular dynamics software developed at RIKEN, Japan. It is designed for large-scale biomolecular simulations on supercomputers, featuring excellent parallel scaling and support for enhanced sampling methods.

**Scientific domain**: Large-scale biomolecular simulations, supercomputing  
**Target user community**: Researchers running large-scale MD on supercomputers

## Theoretical Methods
- Classical molecular dynamics
- Replica exchange methods (REMD, gREST)
- Gaussian accelerated MD (GaMD)
- String method
- Multiple force fields (AMBER, CHARMM, GROMOS)
- Coarse-grained models

## Capabilities (CRITICAL)
- Massively parallel MD
- Replica exchange MD
- Gaussian accelerated MD
- String method for pathways
- Coarse-grained simulations
- Cellular-scale simulations
- GPU acceleration

## Key Strengths

### Parallel Performance:
- Excellent weak scaling
- Fugaku supercomputer optimized
- Millions of atoms
- Hybrid MPI/OpenMP

### Enhanced Sampling:
- Multiple REMD variants
- GaMD implementation
- String method
- gREST method

## Inputs & Outputs
- **Input formats**:
  - PDB structures
  - AMBER prmtop
  - CHARMM PSF
  - GROMACS top
  
- **Output data types**:
  - DCD trajectories
  - Restart files
  - Energy logs

## Interfaces & Ecosystem
- **CHARMM-GUI**: System setup
- **AMBER/CHARMM**: Force fields
- **Fugaku**: Optimized for

## Advanced Features
- **gREST**: Generalized REST
- **GaMD**: Gaussian accelerated MD
- **REUS**: Replica exchange umbrella sampling
- **String method**: Reaction pathways
- **Cellular MD**: Large-scale simulations
- **Hybrid parallel**: MPI + OpenMP + GPU

## Performance Characteristics
- Excellent parallel scaling
- Optimized for ARM (Fugaku)
- GPU support
- Efficient for large systems

## Computational Cost
- Scales to 100,000+ cores
- Efficient for large biomolecules
- GPU acceleration available
- Overall: Excellent for supercomputers

## Best Practices
- Use hybrid parallelization
- Choose appropriate REMD method
- Validate with smaller systems first
- Use CHARMM-GUI for setup

## Limitations & Known Constraints
- Supercomputer focus
- Less desktop-friendly
- Steeper learning curve
- Japanese documentation primary

## Application Areas
- Large biomolecular systems
- Virus capsids
- Membrane proteins
- Cellular-scale simulations
- Drug discovery

## Community and Support
- RIKEN development team
- Documentation
- Workshops (Japan)
- GitHub issues

## Verification & Sources
**Primary sources**:
1. Website: https://www.r-ccs.riken.jp/labs/cbrt/
2. J. Jung et al., J. Phys. Chem. B 128, 5028 (2024)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, LGPL)
- Well-documented
