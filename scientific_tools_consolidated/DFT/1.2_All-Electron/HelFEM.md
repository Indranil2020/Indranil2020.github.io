# HelFEM

## Official Resources
- Source Repository: https://github.com/susilehtola/HelFEM
- Documentation: https://github.com/susilehtola/HelFEM/blob/master/README.md
- License: GNU General Public License v3.0

## Overview
HelFEM (Helsinki Finite Element Method) is a software package for performing fully numerical electronic structure calculations on atoms and diatomic molecules. By using the Finite Element Method (FEM), it avoids the bias of basis set limits (like Gaussians or Slater orbitals), providing benchmark-quality all-electron results.

**Scientific domain**: Atomic physics, diatomic molecules, benchmark data generation, numerical method development
**Target user community**: Developers of functionals, researchers needing basis-set-limit reference data

## Theoretical Methods
- Density Functional Theory (DFT)
- Hartree-Fock (HF)
- Finite Element Method (FEM) discretization
- All-Electron full-potential formulation
- Radial grid solvers for atoms
- 2D FEM solvers for diatomics

## Capabilities
- Atomic calculations (spherically symmetric)
- Diatomic molecular calculations
- Hundreds of exchange-correlation functionals (via LibXC)
- Precise total energies and eigenvalues
- Generation of atomic potentials

## Key Strengths
### Numerical Precision
- Free from basis set superposition error (BSSE)
- Converges systematically to the complete basis set limit
- Ideal for benchmarking other codes

### Versatility
- Supports a vast array of functionals
- Flexible FEM grids

## Inputs & Outputs
- **Input**: Command-line arguments and input scripts
- **Output**: High-precision energies, orbital eigenvalues, potential files

## Interfaces & Ecosystem
- **Core Language**: C++ (uses Armadillo linear algebra).
- **Input**: Command-line interface with scriptable input files.
- **Libraries**: Links against LibXC for functionals.

## Advanced Features
- **Basis Set Limit**:
  - Provides reference values free from basis set truncation error.
- **Radial Solvers**:
  - Extremely high precision for atomic systems.
- **Diatomic FEM**:
  - specialized 2D finite element solver for diatomic molecules.

## Performance Characteristics
- **Focus**: Precision over speed.
- **Parallelization**: Standard LAPACK/BLAS threading; not designed for massive scaling like ErgoSCF.

## Community and Support
- **Development**: Active single-developer/small-team project.
- **Updates**: Regular commits on GitHub.

## Computational Cost
- **Scaling**: Depends on grid density; generally more expensive than Gaussian codes but for much higher precision on small systems.
- **Focus**: Accuracy over large-scale system size.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/susilehtola/HelFEM
2. S. Lehtola, "HelFEM: A finite element code for Hartree-Fock and density functional theory calculations on atoms and diatomic molecules"

**Confidence**: VERIFIED
**Status**: Open Source, Research Code
