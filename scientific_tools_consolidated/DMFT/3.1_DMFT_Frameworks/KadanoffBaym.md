# KadanoffBaym.jl

## Official Resources
- Homepage: https://github.com/NonequilibriumDynamics/KadanoffBaym.jl
- Documentation: https://nonequilibriumdynamics.github.io/KadanoffBaym.jl/
- Source Repository: https://github.com/NonequilibriumDynamics/KadanoffBaym.jl
- License: MIT License

## Overview
KadanoffBaym.jl is a Julia language package designed for solving the Kadanoff-Baym equations (KBE) for non-equilibrium Green's functions (NEGF). It provides efficient, adaptive time-stepping solvers for time-dependent many-body problems, allowing researchers to simulate the dynamics of interacting quantum systems out of equilibrium.

**Scientific domain**: Non-equilibrium quantum dynamics, Many-body physics, Ultrafast spectroscopy
**Target user community**: Theorists studying time-evolution of quantum systems, NEGF practitioners

## Theoretical Methods
- Kadanoff-Baym Equations (KBE)
- Non-Equilibrium Green's Functions (NEGF)
- Keldysh formalism (Time-ordered, Anti-time-ordered, Lesser, Greater components)
- Adaptive time-stepping integrators
- 2-time integration schemes
- Mean-field and correlated approximations (e.g., GKBA, Second Born, GW, T-matrix)

## Capabilities (CRITICAL)
- **Time Evolution**: Solves KBE for the full two-time Green's function $G(t, t')$.
- **Adaptive Solvers**: Implements adaptive stepping to handle rapid dynamics efficiently, crucial for initial transients.
- **Symmetries**: Utilization of time-reversal and particle-hole symmetries to reduce storage and compute.
- **Approximations**: Supports various self-energy approximations (HF, Second Born, GW, T-matrix).
- **Observables**: Calculation of time-dependent currents, densities, and spectral properties.

## Key Features

### Julia Implementation:
- High performance akin to C/Fortran.
- Easy extensibility and integration with other Julia packages.
- Generic type support (e.g., for different precision or matrix types).

### Adaptive Integration:
- Uses advanced algorithms to adapt the time step based on the dynamics.
- Primary integrator: `kbsolve!`.

### Keldysh Components:
- Explicit handling of the contour-ordered Green's function components required for non-equilibrium physics.

## Inputs & Outputs
- **Input formats**:
  - Julia scripts defining Hamiltonians and Self-energies.
  - Initial state definitions.
- **Output data types**:
  - `GreenFunction` objects containing $G^<(t,t')$, $G^>(t,t')$, etc.
  - Time traces of observables.

## Interfaces & Ecosystem
- **Julia Physics Ecosystem**: Part of the `NonequilibriumDynamics` organization (Github).
- **Performance**: Can utilize Julia's threading and GPU capabilities (via array types).

## Workflow and Usage
1.  Define the initial Hamiltonian and the interaction self-energy in Julia.
2.  Define the initial Green's function state `u0`.
3.  Call `kbsolve!(fv!, fd!, u0, (t0, tmax))` to propagate the Green's function in time.
    *   `fv!`: In-place function for the volatile part of the equation.
    *   `fd!`: In-place function for the memory-dependent part (self-energy convolution).
4.  Analyze the output `GreenFunction` object, e.g., using Wigner transforms.

## Performance Characteristics
- **Memory**: $O(N_t^2)$ storage for full two-time Green's functions (challenge for long times).
- **Speed**: Optimized Julia code, effective use of BLAS/LAPACK.


## Comparisons with Other Codes
- **vs Jiezi/NEMO5**: These are typically device simulators (NEGF+Poisson); KadanoffBaym.jl focuses on fundamental KBE time-evolution and many-body approximations.
- **vs NESSY**: Both handle time dynamics, but KadanoffBaym.jl leverages Julia's ecosystem and adaptive solvers.
- **vs Static NEGF**: KadanoffBaym.jl is explicitly for *time-dependent* non-equilibrium problems (transients, pulses).

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/NonequilibriumDynamics/KadanoffBaym.jl
2. Documentation: https://nonequilibriumdynamics.github.io/KadanoffBaym.jl/
3. Related Publications: Papers utilizing KadanoffBaym.jl for NEGF simulations.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Language: Native Julia
- Focus: Specialized solver for KBE/NEGF
