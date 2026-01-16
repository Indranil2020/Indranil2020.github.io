# M-SPARC

## Official Resources
- **Homepage**: https://github.com/SPARC-X/M-SPARC
- **Source Repository**: https://github.com/SPARC-X/M-SPARC
- **Developer**: SPARC-X Group (Georgia Tech)
- **License**: GPLv3

## Overview
**M-SPARC** (Matlab-Simulation Package for Ab-initio Real-space Calculations) is a real-space Finite-Difference Density Functional Theory (DFT) code written in MATLAB/Octave. It is designed for rapid prototyping and educational purposes, allowing researchers to study small-to-moderate sized systems without the complexity of compiled HPC codes. It supports both isolated systems (molecules) and periodic systems.

**Scientific domain**: Real-Space DFT, Finite-Differences, Algorithm Development.
**Target user community**: Developers of DFT algorithms, educators, and students.

## Theoretical Methods
- **Method**: Real-Space Finite-Difference (RSFD).
- **Basis**: None (Grid-based).
- **Pseudopotentials**: ONCV (Optimized Norm-Conserving Vanderbilt).
- **Solvers**: Chebychev filtering for SCF.

## Capabilities
- **Systems**: Isolated interactions, Periodic crystals, and Surfaces.
- **Properties**: Ground state energy, forces, stresses.
- **Relativity**: Spin-polarized calculations.

## Key Strengths
- **Simplicity**: MATLAB codebase is highly readable and modifiable.
- **No Basis Set**: Real-space grid avoids basis set superposition errors (BSSE).
- **Prototyping**: Ideal testbed for new real-space algorithms before implementing in C/C++ (SPARC).

## Comparison with Other Codes
- **vs SPARC**: M-SPARC is the MATLAB prototype; SPARC is the high-performance C++ version for HPC.
- **vs PARSEC**: Both are real-space codes; M-SPARC is easier to modify (MATLAB) but slower than PARSEC (Fortran).
## Performance Characteristics
- **Speed**: Significantly slower than the C++ version (SPARC) but sufficient for small-system prototyping.
- **Parallelization**: Limited by MATLAB's parallel computing (e.g., `parfor`), less efficient than MPI at scale.

## Limitations & Known Constraints
- **Dependencies**: Requires a MATLAB license (or Octave, though compatibility varies).
- **Scale**: Not intended for large-scale production runs (use SPARC for that).

## Best Practices
- **Prototyping**: implementing a new functional or solver? Do it here first in 100 lines of MATLAB before writing 1000 lines of C++.
- **Validation**: Use M-SPARC results to validate the correctness of SPARC runs on supercomputers.

## Community and Support
- **Support**: Maintained by the SPARC-X group at Georgia Tech.
- **Documentation**: Manuals available within the repository.
## Verification & Sources
**Primary sources**:
1.  **Repository**: [M-SPARC GitHub](https://github.com/SPARC-X/M-SPARC)
2.  **Literature**: Xu, Q., et al. "M-SPARC: A MATLAB-based package..." *SoftwareX* (2020).

**Verification status**: ✅ VERIFIED
