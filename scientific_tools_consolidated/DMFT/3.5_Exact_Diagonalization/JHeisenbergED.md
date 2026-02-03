# JHeisenbergED

## Official Resources
- **Repository**: https://github.com/RudSmo/JHeisenbergED
- **License**: MIT License

## Overview
**JHeisenbergED** is a Julia module designed for the **Exact Diagonalization (ED)** and **time-evolution** of the 1D quantum Heisenberg model. It provides a lightweight, pure-Julia implementation for constructing Hamiltonians of spin-1/2 chains and studying their static (ground state) and dynamic properties. It is particularly useful for pedagogical purposes and for researchers needing a quick, programmable solver for small spin systems.

**Scientific domain**: Quantum Magnetism, Many-Body Dynamics
**Target user community**: Students and Julia users studying spin chains

## Theoretical Methods
- **Exact Diagonalization**: Full diagonalization of the Hamiltonian matrix to obtain all eigenvalues and eigenvectors.
- **Model**: 1D Spin-1/2 Heisenberg Chain ($H = J \sum \mathbf{S}_i \cdot \mathbf{S}_{i+1}$).
- **Time Evolution**: Calculation of $|\psi(t)\rangle = e^{-iHt} |\psi(0)\rangle$ using the full unitary operator.

## Capabilities
- **Hamiltonian Construction**: Efficient generation of sparse matrices for arbitrary chain lengths $L$.
- **Spectral Analysis**: Ground state energy and full spectrum.
- ** dynamics**: Exact simulation of quantum quenches and spin dynamics.
- **Observables**: Magnetization, spin correlations.

## Key Strengths
- **Julia Native**: Leverages Julia's speed and easy syntax, avoiding the "two-language problem" (C++ backend + Python frontend) common in other ED codes.
- **Simplicity**: Minimal codebase focuses on doing one thing (Heisenberg model) well.
- **Dynamics**: Built-in support for time evolution, which is often an advanced feature in larger libraries.

## Inputs & Outputs
- **Inputs**:
  - Chain length $L$.
  - Coupling $J$.
  - Initial state vector.
- **Outputs**:
  - Energies.
  - Time-dependent wavefunctions.

## Interfaces & Ecosystem
- **Dependencies**: Julia `LinearAlgebra`, `SparseArrays`.
- **Integration**: Can be easily combined with other Julia packages (e.g., `Plots.jl`, `DifferentialEquations.jl`).

## Performance Characteristics
- **Scaling**: Limited by exponential Hilbert space ($2^L$). Practical limit $L \approx 14-16$ for full diagonalization on a desktop.
- **Efficiency**: Julia implementation compares favorably to C++ for these sizes.

## Comparison with Other Codes
- **vs. QuSpin**: QuSpin is a mature Python ecosystem handling many models and symmetries. JHeisenbergED is a specific tool for Heisenberg chains in Julia.
- **vs. EDLib**: EDLib is a C++ library; JHeisenbergED is for Julia users.

## Application Areas
- **Quantum Quenches**: Studying thermalization and equilibration in closed quantum systems.
- **Education**: demonstrating ED concepts in a high-level language.

## Community and Support
- **Development**: RudSmo (GitHub).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/RudSmo/JHeisenbergED](https://github.com/RudSmo/JHeisenbergED)
- **Verification status**: ✅ VERIFIED
  - Functional Julia package.
