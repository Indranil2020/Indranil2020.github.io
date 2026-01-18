# MonteCarlo.jl

## Official Resources
- Homepage: https://github.com/carstenbauer/MonteCarlo.jl
- Documentation: https://carstenbauer.github.io/MonteCarlo.jl/stable/
- Source Repository: https://github.com/carstenbauer/MonteCarlo.jl
- License: MIT License

## Overview
**MonteCarlo.jl** is a versatile and extensible Julia library for performing classical and quantum Monte Carlo simulations. It provides a unified framework for simulating a wide variety of physical models, including spin systems (Ising, Heisenberg) and interacting fermions (Hubbard model). Designed with modularity in mind, it allows users to easily define custom lattices, Hamiltonians, and observables while leveraging Julia's performance.

**Scientific domain**: Statistical Physics, Condensed Matter, Quantum Many-Body Physics
**Target user community**: Researchers and students modeling spin/fermion systems on lattices

## Theoretical Methods
- **Classical MC**: Metropolis-Hastings for classical spins (Ising, Potts, Heisenberg).
- **Quantum MC**: Stochastic Series Expansion (SSE) for quantum spins/bosons.
- **Determinant QMC (DQMC)**: Auxiliary-field QMC for fermionic Hubbard models.
- **Lattice Theory**: Supports arbitrary lattice graphs (Chain, Square, Cubic, Honeycomb, etc.).

## Capabilities (CRITICAL)
- **Multi-Model Support**: Integrated solvers for Classical Spins, Quantum Spins, and Fermions.
- **Customizability**: deeply modular architecture allows defining new models effectively.
- **Observables**: Built-in measurements for magnetization, energy, specific heat, Green's functions, and correlation functions.
- **Stabilization**: DQMC implementation uses `safe_mult` parameter for controlling numerical stabilization frequency (UDT decomposition).
- **Checkerboard Decomposition**: Optimized update schemes for short-range interactions.

## Key Features

### Unified Interface:
- Describe Model, Lattice, and Algorithm independently.
- consistent API for setting up and running simulations across different physics domains.

### Julia Performance:
- Pure Julia implementation allowing compiler optimizations across user-defined models.
- Comparison with ALPS and other C++ codes shows competitive performance.

### Parallelism:
- Native Julia parallel computing (Binning analysis, parallel tempering).

## Inputs & Outputs
- **Input formats**:
  - Julia scripts creating `Model`, `Lattice`, and `Simulation` objects.
- **Output data types**:
  - JLD2/HDF5 files storing measured observables and statistics.
  - On-the-fly binary output.

## Interfaces & Ecosystem
- **Lattices**: Can generate standard lattices or import custom ones.
- **Analysis**: Tools for binning analysis and error estimation included.

## Workflow and Usage
```julia
using MonteCarlo
model = HubbardModel(L=4, U=4.0)
sim = Simulation(model, algorithm=DQMC())
run!(sim)
```
The workflow emphasizes scriptability and interactivity within the Julia environment.

## Performance Characteristics
- **Speed**: Optimized linear algebra for DQMC; efficient memory access for spin updates.
- **Flexibility**: Trade-off balanced towards extensibility not just raw speed for one specific model.

## Comparison with Other Codes
| Feature | MonteCarlo.jl | ALPS | ALF |
| :--- | :--- | :--- | :--- |
| **Language** | Julia | C++ | Fortran |
| **Focus** | Multi-method (Classical/Quantum) | Broad library of algorithms | Lattice Fermions (DQMC) |
| **Extensibility** | High (Julia multiple dispatch) | Moderate (C++ templates) | Moderate (Module based) |
| **Performance** | Competitive | High | High |
| **User Base** | Growing (Julia community) | Large established | Specialized |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/carstenbauer/MonteCarlo.jl
2. Documentation: https://carstenbauer.github.io/MonteCarlo.jl/

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Activity: Contributors led by C. Bauer.
- Scope: Broad multi-method MC library.
