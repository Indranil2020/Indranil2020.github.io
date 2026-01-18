# libNEGF

## Official Resources
- **Repository**: https://github.com/libnegf/libnegf
- **License**: LGPL-3.0

## Overview
**libNEGF** is a high-performance **Fortran 2008** library capable of solving the Non-Equilibrium Green's Function (NEGF) equations for quantum transport. It is designed as a modular solver (middleware) that can be embedded into various electronic structure codes (like **OMEN**, **CP2K**, or custom drivers). libNEGF provides efficient algorithms for computing the Green's functions, self-energies, and transmission matrices for both orthogonal and non-orthogonal basis sets.

**Scientific domain**: Quantum Transport, Electronic Structure, High-Performance Computing
**Target user community**: Developers of transport codes and researchers needing a robust NEGF kernel

## Theoretical Methods
- **NEGF Formalism**: Stationary (steady-state) transport calculations.
- **Recursive Green's Functions (RGF)**: specialized algorithms (including block-tridiagonal inversion) for wire-like systems to achieve linear scaling $O(N)$.
- **Sancho-Rubio Method**: Iterative calculation of surface Green's functions for semi-infinite leads.
- **Scattering**: Frameworks for including electron-phonon and electron-photon scattering self-energies (approximated).

## Capabilities
- **Observables**:
  - Transmission Probability $T(E)$.
  - Density of States (DOS).
  - Current density.
- **System Types**:
  - Two-terminal devices (Source-Channel-Drain).
  - Nanowires, FinFETs, 2D materials (within a transport setup).
- **Basis Agnostic**: Works with block-sparse matrices provided by the caller (DFT or Tight-Binding).

## Key Strengths
- **Performance**: Heavily optimized for modern supercomputers, utilizing MPI and OpenMP hybrid parallelism.
- **Scalability**: Can handle systems with >100,000 atoms when using RGF.
- **Modularity**: Clean API separates the physics of the Hamiltonian (handled by the caller) from the matrix algebra of transport (handled by libNEGF).

## Inputs & Outputs
- **Inputs**:
  - Hamiltonians ($H$) and Overlap ($S$) matrices (CSR format).
  - Energy grid.
  - Contact self-energies (or descriptions to compute them).
- **Outputs**:
  - Green's function matrices ($G^R, G^<, G^>$).
  - Transmission scalars.
  - Parallel distributed data structures.

## Interfaces & Ecosystem
- **Integrated into**:
  - **OMEN**: A state-of-the-art atomistic quantum transport solver.
  - **CP2K**: Used for transport calculations in the CP2K suite.
- **Usage**: Typically called from Fortran or C/C++ codes.

## Performance Characteristics
- **Speed**: Optimized linear algebra using BLAS/LAPACK and ScaLAPACK.
- **Parallelism**: Multiple levels of parallelism (Energy points, Spatial domain decomposition, Linear algebra blocks).

## Comparison with Other Codes
- **vs. TB_Sim (TBT)**: Similar low-level solver role, but libNEGF is open-source and modern Fortran.
- **vs. Kwant**: Kwant is a user-facing Python package for models; libNEGF is a backend library for massive atomistic simulations.

## Community and Support
- **Development**: Developed at ETH Zurich (Integrated Systems Laboratory) and partners.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/libnegf/libnegf](https://github.com/libnegf/libnegf)
- **Primary Publication**: Check OMEN/libNEGF citations from ETH researchers (Mathieu Luisier et al.).
- **Verification status**: ✅ VERIFIED
  - Active and foundational library for HPC transport.
