# XDiag (formerly exactdiag)

## Official Resources
- **Homepage**: https://awietek.github.io/xdiag/
- **Source Repository**: https://github.com/awietek/xdiag
- **Documentation**: https://awietek.github.io/xdiag/
- **Julia Package**: https://github.com/awietek/XDiag.jl
- **License**: Apache License 2.0

## Overview
XDiag (previously known as exactdiag) is a modern, high-performance C++ library for the Exact Diagonalization (ED) of quantum many-body systems. It is designed to be user-friendly while maintaining high efficiency. It provides a Julia interface (XDiag.jl) which makes it accessible for rapid prototyping and scripting. It specializes in spin systems, t-J models, and Hubbard models.

**Scientific domain**: Quantum Many-Body Physics, Exact Diagonalization  
**Target user community**: Theoretical Physicists, Condensed Matter Researchers

## Capabilities (CRITICAL)
- **Models**: Heisenberg, t-J, Hubbard, and various custom spin/fermion Hamiltonians.
- **Symmetries**: Utilization of translation, point group, and spin symmetries to reduce Hilbert space dimension.
- **Algorithms**:
  - Full diagonalization (Lapack)
  - Iterative diagonalization (Lanczos) for ground state and low-lying excited states.
  - Time evolution.
  - Finite temperature calculations (typicality).
- **Interfaces**: C++17 library and a rich Julia interface.
- **Parallelism**: Shared memory parallelization (OpenMP) and distributed memory (MPI).

## Inputs & Outputs
- **Input formats**: Hamiltonian definition via API (C++ or Julia), lattice geometry.
- **Output data types**: Eigenenergies, eigenstates, correlation functions, spectral functions.

## Performance Characteristics
- Highly optimized for modern architectures.
- Can handle Hilbert spaces up to ~50 billion states (with MPI).
- Efficient implementation of matrix-vector multiplications.

## Application Areas
- Frustrated magnetism.
- Strongly correlated electron systems.
- Quantum phase transitions.
- Benchmarking quantum simulators.

## Verification & Sources
- **Source**: Active GitHub repository (awietek/xdiag).
- **Documentation**: Comprehensive and modern.
- **Identity**: Confirmed renaming from `exactdiag` to `XDiag`.

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Exact Diagonalization
