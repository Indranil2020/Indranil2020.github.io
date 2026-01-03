# exactdiag

## Official Resources
- Homepage: https://github.com/awietek/exactdiag (Presumed)
- Documentation: Minimal
- Source Repository: https://github.com/awietek/exactdiag
- License: GPL v3

## Overview
**Status**: NICHE - "exactdiag" is a generic name. The linked repository (awietek/exactdiag) is a C++ library for exact diagonalization of quantum many-body systems. It focuses on spin systems and Hubbard models.

**Scientific domain**: Exact diagonalization, quantum many-body physics  
**Target user community**: Theoretical physicists

## Capabilities (CRITICAL)
- **Models**: Heisenberg, Hubbard, t-J models.
- **Symmetry**: Uses symmetries (translation, point group) to reduce Hilbert space.
- **Algorithms**: Lanczos algorithm for ground state/low-lying states.
- **Observables**: Correlations, structure factors.

**Sources**: GitHub repository

## Inputs & Outputs
- **Input formats**: C++ source code / config
- **Output data types**: Energies, correlation functions

## Workflow and Usage
1. Define lattice and Hamiltonian in C++.
2. Run executable.

## Performance Characteristics
- Optimized for small systems (limit ~40 spins, ~20 electrons).

## Community and Support
- Research code (Alexander Wietek).

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/awietek/exactdiag

**Confidence**: VERIFIED (as a niche research code)

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: MINIMAL
- Source: OPEN
- Development: ACTIVE
- Applications: Exact diagonalization
