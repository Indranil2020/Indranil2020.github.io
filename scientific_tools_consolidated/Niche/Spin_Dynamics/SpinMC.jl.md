# SpinMC.jl

## Official Resources
- **Homepage**: https://github.com/fbuessen/SpinMC.jl
- **Repository**: https://github.com/fbuessen/SpinMC.jl
- **License**: MIT License

## Overview
**SpinMC.jl** is a Julia package dedicated to **classical Monte Carlo** simulations of lattice spin models. It offers a straightforward interface for defining custom unit cells and interaction matrices, enabling the study of Heisenberg, XY, Ising, Dzyaloshinskii-Moriya, and Kitaev interactions on arbitrary lattices. It is designed to efficiently calculate thermodynamic properties and identify magnetic phase transitions.

**Scientific domain**: Classical Statistical Mechanics, Magnetism
**Target user community**: Researchers studying phase transitions in frustrated magnets

## Theoretical Methods
- **Classical Monte Carlo**: Sampling of the partition function $Z = \int d\mathbf{S} e^{-\beta H(\mathbf{S})}$.
- **Algorithms**:
  - **Metropolis-Hastings**: Local spin updates.
  - **Parallel Tempering (Replica Exchange)**: Efficient sampling of systems with complex energy landscapes (e.g., spin glasses, frustrated systems) by swapping configurations between different temperatures.
- **Interactions**: $H = \sum_{ij} \mathbf{S}_i \mathbf{J}_{ij} \mathbf{S}_j - \sum_i \mathbf{h} \cdot \mathbf{S}_i$.

## Capabilities
- **Model Definition**:
  - Arbitrary crystal lattices.
  - Full $3 \times 3$ interaction matrices (generalized exchange).
  - Single-ion anisotropy.
- **Simulation**:
  - Thermalization sweeps.
  - Measurement sweeps.
  - MPI-parallelized Replica Exchange.
- **Observables**:
  - Internal Energy $E$, Specific Heat $C_v$.
  - Magnetization $M$, Susceptibility $\chi$.
  - Static Structure Factor $S(\mathbf{q})$.
  - Binder Cumulants (for criticality).

## Key Strengths
- **Simplicity**: Lower barrier to entry than larger frameworks; models are defined in pure Julia scripts.
- **Frustration Friendly**: Parallel tempering is built-in, which is essential for getting correct results in frustrated systems like the Kagome antiferromagnet.
- **Julia Scalability**: Efficiently scales from a laptop (single core) to a cluster (MPI) without changing code.

## Inputs & Outputs
- **Inputs**: Julia script with `UnitCell`, `Lattice`, and interaction definitions.
- **Outputs**:
  - HDF5 files with time-series data.
  - Post-processed binary files with mean/error of observables.

## Interfaces & Ecosystem
- **Dependencies**: `MPI.jl`, `HDF5.jl`.
- **Ecosystem**: Can interface with plotting libraries for structure factors.

## Performance Characteristics
- **Speed**: Comparable to C++ codes for local updates due to Julia's LLVM compilation.
- **Scalability**: MPI Parallel Tempering enables simulations of large lattices ($L \sim 100$) near critical points.

## Comparison with Other Codes
- **vs. ALPS/MC**: ALPS is the classic standard. SpinMC.jl is a modern, lightweight alternative that is easier to modify and integrates better with Julia Workflows.
- **vs. Sunny.jl**: Sunny.jl is more general (SU(N) states, dynamics). SpinMC.jl is specifically optimized for *classical* thermodynamics of dipoles.

## Application Areas
- **Spin Liquids**: Searching for lack of order in Kitaev/Heisenberg models.
- **Critical Exponents**: finite-size scaling analysis of phase transitions.

## Community and Support
- **Development**: F. Buessen (University of Toronto / KAIST).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/fbuessen/SpinMC.jl](https://github.com/fbuessen/SpinMC.jl)
- **Verification status**: ✅ VERIFIED
  - Active research code.
