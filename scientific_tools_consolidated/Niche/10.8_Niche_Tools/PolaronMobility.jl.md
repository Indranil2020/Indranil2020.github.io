# PolaronMobility.jl

## Official Resources
- **Repository**: https://github.com/Frost-group/PolaronMobility.jl
- **Documentation**: https://github.com/Frost-group/PolaronMobility.jl (README)
- **License**: MIT License

## Overview
**PolaronMobility.jl** is a Julia package dedicated to calculating polaron properties—specifically mobility and effective mass—in polar semiconductors and ionic crystals. It implements Feynman's variational path-integral approach to the **Fröhlich polaron** problem, which describes the interaction of an electron with macroscopic optical phonons. It also extends to the **Holstein polaron** model, making it a versatile tool for studying charge transport limits in materials like halide perovskites and organic semiconductors.

**Scientific domain**: Carrier Transport, Polarons, Condensed Matter
**Target user community**: Researchers studying mobility limits in polar materials

## Theoretical Methods
- **Feynman Path Integral**: Variational solution to the Fröhlich Hamiltonian.
- **Finite Temperature**: Calculation of free energies and thermal averages at $T > 0$.
- **Mobility Theories**:
  - **FHIP**: Feynman-Hellwarth-Iddings-Platzman (limit of low temperature/weak coupling).
  - **Kadanoff**: Boltzmann equation approach suitable for intermediate temperatures.
- **Holstein Model**: Small polaron physics (lattice discreteness).

## Capabilities
- **Observables**:
  - DC Mobility $\mu(T)$.
  - Optical Conductivity $\sigma(\omega)$ (AC response).
  - Effective Mass $m^*$.
- **Materials**: Input of $\epsilon_0, \epsilon_\infty$, and $\omega_{LO}$ allows simulation of specific compounds (e.g., CsPbBr$_3$).
- **Analysis**:
  - Crossover from large to small polarons.

## Key Strengths
- **Speed**: Variational minimization is computationally inexpensive compared to Monte Carlo or Green's function methods, allowing near-instant results for material screening.
- **Material-Specific**: Designed to take material parameters directly, bridging model physics with real material contexts.
- **Julia Implementation**: Clean code that is easy to integrate into larger material screening pipelines.

## Inputs & Outputs
- **Inputs**: Dielectric constants, Phonon frequency, Effective mass (band mass).
- **Outputs**: Mobility values, plots of $\mu$ vs $T$.

## Interfaces & Ecosystem
- **Dependencies**: Julia scientific stack (`Optim.jl`, `QuadGK.jl`).

## Performance Characteristics
- **Efficiency**: Very high. Solving the variational equations takes milliseconds.
- **Scalability**: trivially parallelizable over different materials or temperatures.

## Comparison with Other Codes
- **vs. EPW**: EPW is a first-principles code that calculates $g_{mn\nu}(k,q)$. PolaronMobility.jl uses a continuum model (Fröhlich). EPW is more accurate for specific band structures; PolaronMobility.jl captures the non-perturbative polaron state better in the continuum limit.
- **vs. DiagMC**: Diagrammatic Monte Carlo is exact but slow. This code is approximate (variational) but fast.

## Application Areas
- **Halide Perovskites**: Explaining the modest mobility and "phonon glass" behavior.
- **Organics**: Understanding transport in soft polar lattices.

## Community and Support
- **Development**: Frost Group (Imperial College London).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/Frost-group/PolaronMobility.jl](https://github.com/Frost-group/PolaronMobility.jl)
- **Verification status**: ✅ VERIFIED
  - Active research tool.
