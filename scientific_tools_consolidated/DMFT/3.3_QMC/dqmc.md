# dqmc (Carsten Bauer)

## Official Resources
- Homepage: https://github.com/carstenbauer/dqmc
- Documentation: https://github.com/carstenbauer/dqmc (README/Examples)
- Source Repository: https://github.com/carstenbauer/dqmc
- License: MIT License

## Overview
**dqmc** is a high-performance Julia-based Determinant Quantum Monte Carlo (DQMC) code developed by Carsten Bauer. It is specifically tailored for simulating strongly correlated electron systems, particularly focusing on quantum critical points in metals (e.g., antiferromagnetic quantum critical point). It leverages the power of Julia and C++ to achieve performance competitive with traditional Fortran/C++ implementations while maintaining code readability and extensibility.

**Scientific domain**: Strongly Correlated Electrons, Quantum Criticality, Lattice Models
**Target user community**: Condensed matter theorists studying 2D lattice fermions and quantum phase transitions

## Theoretical Methods
- **Determinant QMC**: Finite-temperature auxiliary-field QMC for fermions on lattices.
- **Hubbard-Stratonovich Transformation**: Decouples interactions to quadratic forms coupled to bosonic fields.
- **Stabilization Algorithms**: Implements robust stabilization schemes (UDT/QR decompositions) to handle the sign problem and numerical instabilities at low temperatures.

## Capabilities (CRITICAL)
- **Models**: Efficient simulation of Hubbard and Spin-Fermion models.
- **Lattices**: Optimized for 2D square lattices, scalable to large system sizes.
- **Numerical Stability**: Advanced stabilization routines (e.g., from `StableDQMC.jl`) to ensure accurate Green's function computation.
- **Parallelization**: Supports MPI and Julia's threading for running parallel Markov chains or parallel updates.

## Key Features

### High Performance:
- **Optimization**: Critical kernels are highly optimized in Julia or call out to C++ backends.
- **Comparable to C++**: Benchmarks show performance matching or exceeding dedicated C++ codes.

### Algorithmic robustness:
- **Pivoted QR**: Uses pivoted QR decompositions (`StableDQMC.jl`) for numerical stabilization, often superior to standard SVD.
- **LOH Method**: Implementing stable matrix inversion (Linear Operator H).

### Julia Integration:
- **Ecosystem**: Part of a suite of tools including `MonteCarlo.jl` and `StableDQMC.jl`.
- **Interactivity**: Can be run and analyzed interactively in Julia REPL or notebooks.

## Inputs & Outputs
- **Input formats**:
  - Julia scripts or parameter files defining model constants (U, t, beta, etc.).
- **Output data types**:
  - Binary/HDF5 files containing Green's functions, susceptibilities, and observables.
  - Logs of Monte Carlo acceptance rates and measurements.

## Interfaces & Ecosystem
- **Dependencies**: `StableDQMC.jl` (linear algebra), MPI.jl.
- **Analysis**: Data structure compatible with Julia's data analysis tools.

## Workflow and Usage
Users typically define a simulation script in Julia specifying the lattice size, temperature, interaction strength, and number of sweeps. The `dqmc` solver is called to perform the thermalization and measurement steps, outputting correlation functions for post-processing.

## Performance Characteristics
- **Speed**: Highly optimized for specific model classes (quantum critical metals).
- **Scalability**: Good scaling with system size ($N^3$) and parallelization over walkers.

## Comparison with Other Codes
| Feature | dqmc (Julia) | ALF | QUEST |
| :--- | :--- | :--- | :--- |
| **Language** | Julia | Fortran | Fortran |
| **Stabilization** | Pivoted QR / UDT | UDT | UDT / CP |
| **Flexibility** | High (Scriptable) | Moderate | Moderate |
| **Performance** | High (near C++) | High | High |
| **Focus** | Quantum Criticality / 2D | General Lattice | Hubbard Ground State |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/carstenbauer/dqmc
2. "StableDQMC.jl" Repository: https://github.com/carstenbauer/StableDQMC.jl
3. Publications by C. Bauer on quantum critical metals using this code.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Language: Julia (with high-performance focus)
- niche: High-precision DQMC for specific 2D problems.
