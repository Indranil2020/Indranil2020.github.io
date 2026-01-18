# xdiag

## Official Resources
- Homepage: https://github.com/awietek/xdiag
- Documentation: https://awietek.github.io/xdiag/
- Source Repository: https://github.com/awietek/xdiag (and https://github.com/awietek/XDiag.jl)
- License: Apache License 2.0

## Overview
**xdiag** is a modern, high-performance software package for the **Exact Diagonalization (ED)** of quantum many-body systems. Created by Alexander Wietek, it features a dual-language architecture: a highly optimized C++ core for computational efficiency and a Julia wrapper (`XDiag.jl`) for a user-friendly, high-level interface. It is designed to solve generic spin, boson, and fermion models on arbitrary geometries, exploiting symmetries to reach the largest possible system sizes.

**Scientific domain**: Quantum Many-Body Physics, Strongly Correlated Systems, Magnetism
**Target user community**: Theorists using ED for spectral functions, ground states, and time evolution

## Theoretical Methods
- **Exact Diagonalization**: Full or iterative (Lanczos/Davidson) diagonalization of sparse Hamiltonians.
- **Symmetries**:
  - **Translation**: Momentum conservation ($k$-blocks).
  - **Point Group**: Irreducible representations of lattice symmetries.
  - **U(1) / SU(2)**: Conservation of particle number and spin projection $S_z$.
- **Sublattice Coding**: Efficient state representation for spin systems.

## Capabilities (CRITICAL)
- **Hilbert Spaces**: Supports Spin-1/2, Spin-S, t-J, Fermion (Hubbard), and Boson models.
- **General Hamiltonians**: Flexible input for arbitrary interaction terms (bond-dependent, long-range).
- **Large Systems**: Uses distributed memory (MPI) and shared memory (OpenMP) parallelism to solve systems beyond trivial sizes (e.g., 40+ spins).
- **Time Evolution**: Real-time and Imaginary-time evolution algorithms.
- **Spectral Functions**: Calculation of $A(k, \omega)$ and dynamical structure factors.

## Key Features

### Architecture:
- **C++ Core (`libxdiag`)**: Implements fast state lookup (Lin tables), matrix-vector multiplication, and parallelization.
- **Julia Interface**: Seamless usage from Julia, allowing easy plotting and post-processing.

### Algorithms:
- **Distributed Hashing**: Efficient parallelization of state space across nodes without full replication.
- **Iterative Solvers**: Bindings to fast eigensolvers (e.g., Primme or custom Lanczos implementations).

## Inputs & Outputs
- **Input formats**:
  - Julia scripts defining the `Coupling`, `Site`, and `Geometry`.
  - Input files (TOML/JSON) if using the C++ executable directly.
- **Output data types**:
  - Eigenvalues and Eigenvectors.
  - HDF5 files for heavy data storage.

## Interfaces & Ecosystem
- **Dependencies**: C++17 compiler, MPI, BLAS/LAPACK.
- **Julia**: `XDiag.jl` registered in Julia General registry.

## Workflow and Usage
**Julia Example**:
```julia
using XDiag
block = Spinhalf(N=16, n_up=8)
ops = OpSum()
ops += "Sz", 1, "Sz", 2
H = Hamiltonian(ops, block)
eigs, vecs = eigsolve(H)
```

## Performance Characteristics
- **Speed**: State-of-the-art implementation using bit-manipulation and optimized lookup tables.
- **Scaling**: Scales to supercomputing clusters via MPI.

## Comparison with Other Codes
| Feature | xdiag | QuSpin | EDKit.jl |
| :--- | :--- | :--- | :--- |
| **Core** | C++ (Distributed) | Python/C++ | Julia |
| **Parallelism** | MPI + OpenMP | OpenMP (Shared) | Threads (Julia) |
| **Scaling** | Cluster / Supercomputer | Single Node | Single Node |
| **Symmetries** | Extensive (Point/Space) | Moderate | U(1) / Translation |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/awietek/xdiag
2. Documentation: https://awietek.github.io/xdiag/
3. Publication: "XDiag: A high-performance exact diagonalization library" (referenced in SciPost/arXiv).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (Apache 2.0)
- Authors: A. Wietek (Flatiron Institute).
- Quality: Professional-grade research software.
