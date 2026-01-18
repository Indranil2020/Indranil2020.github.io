# Sunny.jl

## Official Resources
- **Homepage**: https://github.com/SunnySuite/Sunny.jl
- **Documentation**: https://sunnysuite.github.io/Sunny.jl/dev/
- **Repository**: https://github.com/SunnySuite/Sunny.jl
- **License**: MIT License

## Overview
**Sunny.jl** is a cutting-edge Julia package for the simulation of **spin dynamics** and magnetic properties in crystal systems. What sets it apart is its formulation based on **SU(N) coherent states**, which generalizes the traditional Landau-Lifshitz-Gilbert (LLG) dynamics of dipoles to include multipolar moments (quadrupoles, octupoles) on an equal footing. It unifies Classical Monte Carlo, Molecular Dynamics (LLG), and Linear Spin Wave Theory (LSWT) in a single, high-performance framework.

**Scientific domain**: Quantum Magnetism, Spin Dynamics
**Target user community**: Theorists studying frustrated magnets, neutron scattering, and multipolar orders

## Theoretical Methods
- **SU(N) Coherent States**: Representation of local magnetic moments allowing for on-site anisotropy and multipolar order parameters.
- **Dynamics**: Generalized Landau-Lifshitz-Gilbert (LLG) equation for SU(N) spins.
- **Thermodynamics**: Classical Monte Carlo (Metropolis) and Langevin dynamics.
- **Excitations**: Linear Spin Wave Theory (LSWT) generalized for SU(N) ground states.

## Capabilities
- **Simulations**:
  - Time-evolution of spin textures.
  - Temperature-dependent phase diagrams.
  - Calculation of Magnon dispersion relations $\omega(\mathbf{k})$.
- **Observables**:
  - Static Structure Factor $S(\mathbf{q})$.
  - Dynamical Structure Factor $S(\mathbf{q}, \omega)$ (comparable to Neutron Scattering).
- **Interactions**:
  - Heisenberg, Dzyaloshinskii-Moriya, Single-Ion Anisotropy.
  - Biquadratic and higher-order exchange.
  - Long-range dipole-dipole interactions.

## Key Strengths
- **Beyond Dipoles**: The ability to treat $S > 1/2$ systems using SU(N) logic allows it to capture physics (like quadrupolar ordering) that standard codes (SpinW, UppASD) miss.
- **Performance**: Leveraging Julia, it matches or exceeds the performance of C++ codes while maintaining script-level flexibility.
- **Unified Workflow**: One code for finding the ground state (MC), relaxing it (LLG), and calculating the spectrum (LSWT).

## Inputs & Outputs
- **Inputs**:
  - Crystal structure (CIF import supported).
  - Exchange interactions $J_{ij}$ (symmetry analysis helps identify allowed terms).
- **Outputs**:
  - Structure factors (NumPy-compatible arrays).
  - Dispersion plots.

## Interfaces & Ecosystem
- **Crystals.jl**: For crystal structure handling.
- **GLMakie.jl**: For interactive 3D visualization of spin configurations.

## Performance Characteristics
- **Speed**: Highly optimized memory layout and loop structure.
- **Parallelism**: Multi-threading for calculating $S(\mathbf{q}, \omega)$ over momentum space.

## Comparison with Other Codes
- **vs. SpinW**: SpinW (MATLAB) describes dipoles. Sunny.jl describes SU(N) coherent states. Sunny is free, open-source, and faster for many tasks.
- **vs. UppASD**: UppASD is a powerful Fortran code for Atomistic Spin Dynamics. Sunny offers a more modern Julia interface and the unique SU(N) feature set.

## Application Areas
- **Neutron Scattering**: Direct simulation of experimental $S(\mathbf{q}, \omega)$.
- **Multipolar Magnetism**: Studying hidden order in heavy fermion / f-electron materials.

## Community and Support
- **Development**: Sunny Suite Team (Johns Hopkins University, various labs).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/SunnySuite/Sunny.jl](https://github.com/SunnySuite/Sunny.jl)
- **Verification status**: ✅ VERIFIED
  - Active, high-profile research code.
