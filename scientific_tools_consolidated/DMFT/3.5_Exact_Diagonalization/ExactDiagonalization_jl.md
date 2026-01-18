# ExactDiagonalization.jl (QuantumLattices)

## Official Resources
- Homepage: https://github.com/Quantum-Many-Body/ExactDiagonalization.jl
- Documentation: https://quantum-many-body.github.io/ExactDiagonalization.jl/dev/
- Source Repository: https://github.com/Quantum-Many-Body/ExactDiagonalization.jl
- License: MIT License

## Overview
**ExactDiagonalization.jl** is a comprehensive Julia package for the exact diagonalization of quantum many-body systems. It serves as the ED engine within the **QuantumLattices.jl** ecosystem. It provides a unified interface to solve bosonic, fermionic, and spin models defined using the symbolic operator algebra of QuantumLattices. It supports both full diagonalization (for small systems/thermal states) and iterative methods (Lanczos/Arnoldi) for ground states and low-lying excitations.

**Scientific domain**: Condensed Matter Theory, Quantum Magnetism, Strongly Correlated Systems
**Target user community**: Users of the QuantumLattices.jl framework

## Theoretical Methods
- **Symbolic Operator Mapping**: Automatic conversion of symbolic Hamiltonians (from QuantumLattices) to matrix representations.
- **Matrix Diagonalization**: Wrapper around Julia's `Eigen` (dense) and `Arpack` (sparse) solvers.
- **Thermodynamics**: Calculation of thermal averages via full spectrum or trace estimations.
- **Dynamical Correlations**: Kubo formula and Lehmann representation for structure factors and susceptibilities.

## Capabilities (CRITICAL)
- **Generic Models**: Supports Spins, Fermions, Bosons, and mixtures.
- **Unified Workflow**: Seamless integration with the lattice generation and model definition of QuantumLattices.jl.
- **Observables**: Automated calculation of local observables and correlation functions.
- **Path Integral**: Some support for path integral formulations (referenced in ecosystem).

## Key Features

### Ecosystem Integration:
- **QuantumLattices.jl**: Relies on this for geometry (UnitCell, Lattice) and Algebra (Hilbert space).
- **Extensibility**: Part of a modular suite including DMRG and QMC codes under `Quantum-Many-Body`.

### Algorithms:
- **ED**: Full diagonalization for thermodynamic properties.
- **Krylov**: Iterative solvers for large sparse matrices.

## Inputs & Outputs
- **Input formats**:
  - Julia scripts using `Lattice`, `Hilbert`, and `Hamiltonian` objects from QuantumLattices.
- **Output data types**:
  - `Eigen` objects (values/vectors).
  - DataFrames or Dictionaries of results.

## Interfaces & Ecosystem
- **Dependencies**: `QuantumLattices.jl`, `SparseArrays`, `LinearAlgebra`, `KrylovKit`.
- **Sister Packages**: `DMRG.jl`, `QuantumMC.jl`.

## Workflow and Usage
```julia
using QuantumLattices, ExactDiagonalization
lattice = Lattice(...)
hilbert = Hilbert(site=>Spin{1//2}() for site in lattice)
H = Hamiltonian(...)
ed = ED(lattice, hilbert, H)
solve!(ed)
```

## Performance Characteristics
- **Speed**: Efficient construction of matrices from symbolic rules; standard solver performance.
- **Versatility**: Prioritizes generality over extreme optimization for one specific model type.

## Comparison with Other Codes
| Feature | ExactDiagonalization.jl | EDKit.jl | KrylovKit.jl |
| :--- | :--- | :--- | :--- |
| **Scope** | Physics Engine (Spins/Fermions) | Physics Engine (General) | Generic Linear Algebra Solver |
| **Integration** | Tight (QuantumLattices.jl) | Standalone | used by ED codes |
| **User Level** | High (Physics Objects) | High (Physics Objects) | Low (Matrices/Vectors) |
| **Ecosystem** | Quantum-Many-Body | Roger-luo/Yao | JuliaLinearAlgebra |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/Quantum-Many-Body/ExactDiagonalization.jl
2. Documentation: https://quantum-many-body.github.io/ExactDiagonalization.jl/dev/

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Ecosystem: Core component of Quantum-Many-Body organization.
