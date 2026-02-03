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

## Best Practices
- Use appropriate symmetry function parameters
- Validate on test set
- Consider AmpTorch for GPU acceleration
- Start with default descriptors

## Limitations & Known Constraints
- CPU-based (AmpTorch for GPU)
- Older architecture (pre-equivariant era)
- Less accurate than equivariant methods
- Symmetry function choice critical

## Application Areas
- Surface science
- Catalysis
- Materials modeling
- Educational purposes

## Comparison with Other Codes
- **vs AmpTorch**: AMP CPU-based original, AmpTorch GPU PyTorch version
- **vs NequIP/MACE**: AMP descriptor-based, others equivariant (more accurate)
- **vs N2P2**: Both Behler-Parrinello style, different implementations
- **Unique strength**: Simplicity, good documentation, educational value, ASE integration

## Community and Support
- Established codebase
- Good documentation
- ASE community
- Published tutorials

## Verification & Sources
**Primary sources**:
1. Documentation: https://amp.readthedocs.io/
2. A. Khorshidi & A.A. Peterson, Comput. Phys. Commun. 207, 310 (2016)

**Secondary sources**:
1. AMP tutorials
2. AmpTorch documentation
3. Published catalysis applications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (Bitbucket, GPL-3.0)
- Academic citations: >500
- Established codebase
