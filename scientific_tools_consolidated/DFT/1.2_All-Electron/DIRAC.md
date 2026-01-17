# DIRAC

## Official Resources
- Homepage: https://www.diracprogram.org/
- Source Repository: https://gitlab.com/dirac/dirac
- Documentation: https://www.diracprogram.org/doc/release-24/
- License: LGPL-2.1

## Overview
DIRAC (Program for Atomic and Molecular Direct Iterative Relativistic All-electron Calculations) is the premier general-purpose relativistic quantum chemistry program. It is designed to treat relativistic effects in molecules with the highest possible accuracy, serving as a benchmark for approximate methods. It enables all-electron calculations using the 4-component Dirac-Coulomb implementation of both Density Functional Theory (DFT) and wavefunction-based correlation methods.

**Scientific domain**: Relativistic quantum chemistry, heavy-element chemistry, hazardous materials (actinides)
**Target user community**: Quantum chemists requiring high-precision relativistic treatments

## Theoretical Methods
- **Hamiltonians**:
    - 4-component Dirac-Coulomb
    - 4-component Dirac-Coulomb-Breit
    - Exact Two-Component (X2C)
    - Levy-Leblond (non-relativistic limit)
- **DFT**: 4-component Dirac-Kohn-Sham (DKS) with a wide range of non-collinear functionals.
- **Wavefunction Methods**: MP2, Coupled Cluster (CCSD(T)), CI, MCSCF at the 4-component level.
- **Basis Sets**: Large library of relativistic basis sets (Dyall, etc.).

## Capabilities
- **Molecular Properties**: Unrivaled accuracy for NMR shieldings, spin-rotation constants, electric field gradients, and parity violation.
- **Excited States**: Relativistic Linear Response TDDFT and EOM-CC.
- **Open-Shell Systems**: Sophisticated treatment of open-shell species including spin-orbit coupling.
- **Solvation**: PCM and explicit solvent models suitable for relativistic calculations.

## Key Strengths
### The Relativistic Benchmark
- DIRAC is often used to validate results from more approximate (2-component or scalar relativistic) codes. Its 4-component treatment is the "gold standard" for molecular relativity.

### Comprehensive Methodology
- Uniquely combines high-level wavefunction theory (CC, CI) with relativistic DFT in a single package.

## Inputs & Outputs
- **Input**: `mol` (geometry) and `inp` (calculation control) files, or Python scripting interface.
- **Output**: Detailed analysis of relativistic wavefunctions, energies, and property tensors.

## Interfaces & Ecosystem
- **Launcher**: `pam` script for easy execution and memory management.
- **Python**: Exposes a Python API for complex workflows.
- **Parallelization**: MPI for distributed memory, OpenMP for shared memory.
- **QCSchema**: Supports modern JSON-based output standards.

## Computational Cost
- **High**: 4-component calculations are significantly more expensive than non-relativistic ones (often 10x-100x).
- **Optimization**: Uses quaternion algebra and symmetry to reduce cost. X2C module offers cheaper alternatives.

## Verification & Sources
**Primary sources**:
1. Official Website: https://www.diracprogram.org/
2. "The DIRAC code for relativistic molecular calculations" (J. Chem. Phys. 152, 204104 (2020))

## Community and Support
- **Forum**: Google Group for users and developers.
- **Workshop**: Annual "Relativistic Quantum Chemistry" schools using DIRAC.
- **Status**: Active development (DIRAC25 released 2025).

**Confidence**: VERIFIED
**Status**: Active, Open Source
