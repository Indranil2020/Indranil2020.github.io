# PEPS (QuantumLiquids)

## Official Resources
- Homepage: https://github.com/QuantumLiquids/PEPS
- Source Repository: https://github.com/QuantumLiquids/PEPS
- License: Open Source (Header-only core)

## Overview
**PEPS** (Projected Entangled Pair States) by the *QuantumLiquids* group is a high-performance C++ library for simulating 2D strongly correlated electron systems. It implements the **PEPS** tensor network ansatz and uses a **Variational Monte Carlo (VMC)** approach to optimize the tensor elements. This combination allows for the rigorous study of 2D fermionic systems (fPEPS) that are challenging for traditional QMC due to the sign problem.

**Scientific domain**: 2D Quantum Lattice Models, Tensor Networks, Strongly Correlated Electrons
**Target user community**: Researchers studying high-Tc superconductivity, topological order, and 2D magnetism

## Theoretical Methods
- **Projected Entangled Pair States (PEPS)**: A 2D generalization of MPS that naturally captures Area Law entanglement.
- **Variational Monte Carlo (VMC)**: stochastic sampling to evaluate the contraction of the norm $\langle\psi|\psi\rangle$ and observables, avoiding the high cost of exact contraction.
- **Stochastic Reconfiguration (SR)**: Optimization method (similar to Natural Gradient) to update variational parameters.
- **Fermionic PEPS (fPEPS)**: Explicit handling of fermion statistics (swap gates).

## Capabilities (CRITICAL)
- **Fermionic Systems**: Native support for fPEPS, crucial for the Hubbard and t-J models.
- **Optimization**: State-of-the-art VMC optimizers including Adam and Stochastic Reconfiguration.
- **Custom Models**: Plugin architecture for defining new Hamiltonians and lattice geometries.
- **Observables**: Built-in measurement toolkit for energy, correlation functions, and structure factors.
- **Header-Only**: Easy integration into other C++ projects.

## Key Features

### Performance:
- **C++20**: modern C++ implementation for efficiency and safety.
- **Parallelism**: MPI support for parallelizing the Monte Carlo sampling steps.
- **Math Libraries**: Links against high-performance BLAS (Intel MKL / OpenBLAS).

### Flexibility:
- **Lattices**: Built-in support for Square and other common 2D lattices.
- **Extensible**: Plugin interfaces for Models and Updaters.

## Inputs & Outputs
- **Input formats**:
  - C++ driver files or configuration scripts.
- **Output data types**:
  - Energy logs, optimized tensor parameters.
  - Measured observables (Spin/Charge correlations).

## Interfaces & Ecosystem
- **Dependencies**: C++20 compiler, CMake, BLAS/LAPACK, MPI.
- **Integration**: Part of the `QuantumLiquids` methods suite (TensorToolkit).

## Workflow and Usage
Users define the model (e.g., Heisenberg or Hubbard) and the PEPS ansatz (bond dimension D) in C++. The VMC optimizer is then run to minimize the energy.
```cpp
// Example conceptual usage
auto model = HubbardModel(Lx, Ly, U);
auto peps = PEPSState(Lx, Ly, D);
VMCOptimizer optimizer(peps, model);
optimizer.run();
```

## Performance Characteristics
- **Cost**: Contraction is \#P-hard, but VMC sampling makes it manageable ($O(D^6)$ or similar depending on update scheme).
- **Scalability**: Excellent scaling with Monte Carlo samples via MPI.

## Comparison with Other Codes
| Feature | PEPS (QuantumLiquids) | TeNeS | ITensor |
| :--- | :--- | :--- | :--- |
| **Method** | VMC Optimization of PEPS | Imaginary Time Evolution | MPS / DMRG |
| **Dimension** | 2D (NxN) | 2D (Infinite) | 1D (Primary) |
| **Contraction** | Monte Carlo Sampling | Deterministic (CTMRG) | Deterministic (Exact) |
| **Fermions** | Native (fPEPS) | Native | Native |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/QuantumLiquids/PEPS
2. Recent arXiv preprints on fPEPS using VMC/SR methods (e.g. from the authors).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- State: Active development.
