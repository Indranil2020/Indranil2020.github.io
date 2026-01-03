# opendf

## Official Resources
- Homepage: https://github.com/CQMP/opendf
- Documentation: https://github.com/CQMP/opendf/blob/master/README.md
- Source Repository: https://github.com/CQMP/opendf
- License: GNU General Public License v2.0

## Overview
opendf is a condensed matter physics code that solves strongly correlated lattice problems (such as the Hubbard model) in finite dimensions using the dual fermion method. It extends DMFT by including non-local correlations through diagrammatic extensions, providing a systematic way to treat spatial correlations beyond local DMFT approximations.

**Scientific domain**: Strongly correlated systems, dual fermion method, diagrammatic extensions of DMFT  
**Target user community**: Researchers studying non-local correlations in strongly correlated materials

## Theoretical Methods
- Dual fermion (DF) method
- Diagrammatic extensions of DMFT
- Non-local correlation effects
- Ladder approximation
- Hubbard model in finite dimensions
- Spatial correlations beyond DMFT
- Vertex function calculations

## Capabilities (CRITICAL)
- Dual fermion calculations
- Non-local correlation effects
- Hubbard model solutions in 2D and 3D
- Systematic improvements beyond DMFT
- Ladder diagrams and higher-order correlations
- Momentum-dependent self-energies
- Vertex functions
- C++ implementation with ALPSCore libraries
- HDF5 data format
- MPI parallelization

**Sources**: Official opendf repository (https://github.com/CQMP/opendf), confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- Configuration files
- Lattice parameters
- Interaction parameters
- DMFT starting solution
- HDF5 input data

**Output data types**:
- Self-energies (momentum-dependent)
- Green's functions
- Vertex functions
- Observables
- HDF5 archives

## Interfaces & Ecosystem
- **ALPSCore**: Built on ALPSCore libraries
- **C++**: Modern C++ implementation
- **HDF5**: Standard data format
- **MPI**: Parallel execution

## Limitations & Known Constraints
- Computational cost higher than DMFT alone
- Requires DMFT solution as starting point
- Complex installation (ALPSCore dependency)
- Expertise in dual fermion method required
- Limited to lattice models
- Documentation technical
- Platform: Linux/Unix

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/CQMP/opendf
2. README and code documentation
3. CQMP group (Correlated Quantum Materials Physics)

**Secondary sources**:
1. Dual fermion method literature
2. Published applications
3. ALPSCore ecosystem
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official repository: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (README)
- Source code: OPEN (GitHub, GPL v2)
- Uses ALPSCore libraries
- Dual fermion method implementation
- Active research code
