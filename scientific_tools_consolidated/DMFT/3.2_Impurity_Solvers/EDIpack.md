# EDIpack

## Official Resources
- Homepage: https://edipack.github.io/EDIpack/
- Documentation: https://edipack.github.io/EDIpack/
- Source Repository: https://github.com/EDIpack/EDIpack
- License: GNU General Public License v3.0

## Overview
EDIpack is a massively parallel exact diagonalization solver for generic quantum impurity problems. It uses Lanczos-based methods to solve multi-orbital impurity models, providing an alternative to quantum Monte Carlo approaches. EDIpack is particularly suited for problems where exact solutions are needed or where QMC suffers from sign problems, and it can handle electron-phonon coupling.

**Scientific domain**: Exact diagonalization, quantum impurity solvers, DMFT  
**Target user community**: Researchers needing exact impurity solutions or studying systems with sign problems

## Theoretical Methods
- Exact diagonalization (ED)
- Lanczos algorithm
- Massively parallel execution
- Multi-orbital Anderson impurity model
- Electron-phonon coupling
- Zero-temperature and finite-temperature
- Superconducting states
- Spin-non-conserving problems

## Capabilities (CRITICAL)
- Exact diagonalization impurity solver
- Multi-orbital quantum impurity problems
- Massively parallel (MPI parallelization)
- Zero-temperature solutions
- Finite-temperature properties
- Electron-phonon coupling support
- Superconducting (s-wave) systems
- Spin-non-conserving Hamiltonians
- Green's functions and observables
- Spectral functions (no analytical continuation needed)
- Real-frequency results
- Integration with DMFT frameworks

**Sources**: Official EDIpack documentation (https://edipack.github.io/EDIpack/), A. Amaricci et al., Comput. Phys. Commun. 273, 108261 (2022)

## Inputs & Outputs
**Input formats**:
- Bath parameters
- Interaction parameters
- Impurity Hamiltonian
- Configuration files

**Output data types**:
- Green's functions (real and imaginary frequencies)
- Self-energies
- Spectral functions
- Observables
- Energy levels and eigenstates
- Partition function

## Interfaces & Ecosystem
- **DMFT frameworks**: Can be used as impurity solver
- **Fortran**: Fortran-based with library interface
- **MPI**: Distributed memory parallelization
- **SciFortran**: Uses SciFortran library

## Limitations & Known Constraints
- Limited to relatively small Hilbert spaces
- Bath discretization required
- Memory intensive for large problems
- Computational cost scales exponentially with system size
- Truncation of Hilbert space necessary for large problems
- Not suitable for very large multi-orbital systems
- Real-frequency only (advantage: no analytical continuation)

## Verification & Sources
**Primary sources**:
1. Official documentation: https://edipack.github.io/EDIpack/
2. GitHub repository: https://github.com/EDIpack/EDIpack
3. A. Amaricci et al., Comput. Phys. Commun. 273, 108261 (2022) - EDIpack paper

**Secondary sources**:
1. EDIpack tutorials and examples
2. Published applications using EDIpack
3. Exact diagonalization literature
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, GPL v3)
- Active development: Regular updates
- Academic citations: Growing user base
- Massively parallel: Significant scalability
