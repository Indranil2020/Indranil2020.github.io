# AMP (Atomistic Machine-learning Package)

## Official Resources
- Homepage: https://amp.readthedocs.io/
- Documentation: https://amp.readthedocs.io/
- Source Repository: https://bitbucket.org/andrewpeterson/amp
- License: GPL-3.0

## Overview
AMP (Atomistic Machine-learning Package) is an open-source package for fitting neural network potentials to atomistic data. It uses Behler-Parrinello symmetry functions as descriptors and provides a simple interface for training potentials from DFT data.

**Scientific domain**: Neural network potentials, atomistic machine learning  
**Target user community**: Researchers training neural network potentials

## Theoretical Methods
- Behler-Parrinello symmetry functions
- Neural network regression
- Gaussian descriptor functions
- Atomic fingerprints

## Capabilities (CRITICAL)
- Neural network potential training
- Symmetry function descriptors
- ASE calculator interface
- Parallel training
- Multiple descriptor types

## Key Strengths

### Simplicity:
- Easy to use
- Good documentation
- ASE integration

### Flexibility:
- Custom descriptors
- Multiple backends

## Inputs & Outputs
- **Input formats**: ASE trajectory, VASP OUTCAR
- **Output data types**: Trained models, energies, forces

## Interfaces & Ecosystem
- **ASE**: Calculator interface
- **AmpTorch**: PyTorch version

## Advanced Features
- **Symmetry functions**: G2, G4 types
- **Parallel training**: MPI support
- **Custom descriptors**: Extensible

## Performance Characteristics
- CPU-based training
- Good for small datasets
- ASE integration efficient

## Computational Cost
- Training: Hours (CPU)
- Inference: Fast
- Overall: Moderate

## Limitations & Known Constraints
- CPU-based (AmpTorch for GPU)
- Older architecture
- Less accurate than equivariant methods

## Application Areas
- Surface science
- Catalysis
- Materials modeling

## Verification & Sources
**Primary sources**:
1. Documentation: https://amp.readthedocs.io/
2. A. Khorshidi & A.A. Peterson, Comput. Phys. Commun. 207, 310 (2016)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (Bitbucket, GPL-3.0)
