# EDKit.jl

## Official Resources
- Homepage: https://github.com/Roger-luo/EDKit.jl
- Source Repository: https://github.com/Roger-luo/EDKit.jl
- License: MIT License

## Overview
**EDKit.jl** is a lightweight Julia package for performing Exact Diagonalization (ED) of generic many-body quantum systems. Developed by Roger Luo (author of Yao.jl), it provides a flexible framework for constructing Hamiltonians, handling symmetries, and performing spectral analysis. It is designed to be extensible, allowing users to define custom operator bases and exploit symmetries like $U(1)$ (particle conservation) and $\mathbb{Z}_2$ (spatial reflection).

**Scientific domain**: Quantum Many-Body Physics, Tensor Networks, Quantum Information
**Target user community**: Researchers prototyping ED codes or studying small quantum lattice systems

## Theoretical Methods
- **Exact Diagonalization**: Construction of sparse Hamiltonian matrices.
- **Symmetry Sectors**: Block-diagonalization using conserved quantities (Total Sz, Particle Number).
- **Matrix-Free Methods**: Compatible with iterative solvers (KrylovKit.jl) for large sparse systems.

## Capabilities (CRITICAL)
- **Hamiltonian Construction**: Efficient generation of Hamiltonians from symbolic operator strings.
- **Symmetry Support**: Built-in support for $U(1)$ symmetry and translation/reflection symmetries (`TranslationalBasis`).
- **Custom Bases**: Extensible architecture to define new local Hilbert spaces and basis mappings.
- **Integration**: Works well with the wider Julia quantum ecosystem (e.g., Yao.jl for Quantum Computing).

## Key Features

### Flexibility:
- **Operator Algebra**: Define models using natural operator language (e.g., `Op("Z", i)`).
- **Generic Types**: Leverages Julia's type system to handle various number types and storage backends.

### Usability:
- **Lightweight**: Minimal dependencies compared to full-blown tensor network frameworks.
- **Pedagogical**: Clean implementation suitable for learning ED techniques.

## Inputs & Outputs
- **Input formats**:
  - Julia scripts defining the lattice, local dimension, and operator terms.
- **Output data types**:
  - Spectrum (Eigenvalues), Eigenstates.
  - Correlation functions.

## Interfaces & Ecosystem
- **Dependencies**: `SparseArrays`, `LinearAlgebra`.
- **Related Tools**: `Yao.jl` (Quantum Computing), `KrylovKit.jl` (Eigensolvers).

## Workflow and Usage
```julia
using EDKit
# Define basis with symmetry
basis = TensorProductBasis(2, 4; qn = TotalSz(0))
# Construct Hamiltonian
H = trans_invariant_hamiltonian(basis, ...)
# Solve
vals, vecs = eigs(H)
```

## Performance Characteristics
- **Speed**: Comparable to C++ for sparse matrix construction due to Julia's efficient loops.
- **Scaling**: Limited by exponential Hilbert space growth (typical ED limit ~20 spins depending on symmetry).

## Comparison with Other Codes
| Feature | EDKit.jl | ExactDiagonalization.jl | xdiag |
| :--- | :--- | :--- | :--- |
| **Language** | Julia | Julia | C++ / Julia Wrapper |
| **Focus** | Lightweight / Prototyping | Integrated (QuantumLattices.jl) | HPC / Production |
| **Symmetry** | U(1), Z2, Translation | General (via QL) | General (Space/Point) |
| **Performance** | Good (Julia) | Good (Julia) | Excellent (C++ Core) |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/Roger-luo/EDKit.jl
2. Author's website/portfolio (Roger Luo).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- State: Stable, utility package.
