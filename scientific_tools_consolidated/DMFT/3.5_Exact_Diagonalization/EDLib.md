# EDLib

## Official Resources
- **Repository**: https://github.com/Q-solvers/EDLib
- **License**: GPL-3.0

## Overview
**EDLib** is a flexible C++ template library for the **Exact Diagonalization (ED)** of quantum many-body systems. It specifically targets fermionic models like the **Hubbard Model** and the **Anderson Impurity Model (AIM)** on finite clusters. Designed with efficiency and modern C++ practices in mind, it provides tools for computing ground state properties, finite-temperature thermodynamics, and spectral functions using Lanczos algorithms.

**Scientific domain**: Strongly Correlated Electrons, DMFT Solvers
**Target user community**: Developers of DMFT codes and theorists studying cluster models

## Theoretical Methods
- **Exact Diagonalization (ED)**: Full or iterative diagonalization of the Hamiltonian matrix.
- **Lanczos Algorithm**:
  - Ground state energy and wavefunction.
  - Finite-Temperature Lanczos (FT-Lanczos) for thermodynamic averages.
- **Green's Functions**:
  - Continued Fraction expansion.
  - Lehmann representation for spectral functions $A(\omega)$.
- **Symmetries**: Utilization of $U(1)$ charge and spin symmetries to block-diagonalize matrices.

## Capabilities
- **Models**:
  - Single and Multi-orbital Hubbard Models.
  - Anderson Impurity Models (AIM) with general bath geometries.
  - t-J Models (via mapping).
- **Observables**:
  - Spectral Functions (DOS).
  - Spin-Spin and Charge-Charge correlations.
  - Local and non-local susceptibilities.
  - Specific Heat and Entropy.

## Key Strengths
- **Library Design**: Header-only style template library makes it easy to include in other larger C++ projects (e.g., as the impurity solver for a DMFT code).
- **Flexibility**: Arbitrary lattice geometries and impurity bath structures can be defined.
- **Performance**: Optimized sparse matrix storage (CSR) and bit-manipulation for fermionic state indexing.

## Inputs & Outputs
- **Inputs**: C++ code defining the model parameters and geometry.
- **Outputs**: Numerical data for spectra and correlation functions.

## Interfaces & Ecosystem
- **Dependencies**: MPI (for cluster parallelism), OpenMP (shared memory), HDF5 (optional for I/O).
- **Integration**: Often used as the backend solver for custom DMFT implementations.

## Performance Characteristics
- **Scaling**: Constrained by the exponential growth of the Hilbert space. feasible for $N \approx 14-16$ sites/orbitals on standard nodes, up to ~20-24 on large clusters with MPI.
- **Parallelism**: Hybrid MPI+OpenMP allows effective utilization of modern HPC nodes.

## Comparison with Other Codes
- **vs. ALPS/ED**: ALPS is a large, integrated application suite. EDLib is a lightweight *library*, offering lower overhead for developers who want to write their own Hamiltonian logic in C++.
- **vs. Hydra**: Hydra is another modern C++ ED library; EDLib has a specific historical focus on impurity models for DMFT.

## Application Areas
- **DMFT**: Solving the effective impurity problem in Dynamical Mean Field Theory.
- **Quantum Dots**: simulating transport and spectra of small interacting dot arrays.
- **Cluster approximations**: DCA (Dynamical Cluster Approximation) studies.

## Community and Support
- **Development**: Q-solvers organization (GitHub).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/Q-solvers/EDLib](https://github.com/Q-solvers/EDLib)
- **Verification status**: ✅ VERIFIED
  - Active modern C++ project.
