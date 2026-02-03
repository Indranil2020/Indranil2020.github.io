# almaBTE

## Official Resources
- **Repository**: https://bitbucket.org/sousaw/almabte
- **Documentation**: http://www.almabte.eu/
- **License**: Apache License 2.0

## Overview
**almaBTE** is a high-performance C++ software package designed for calculating the **lattice thermal conductivity** and other phonon transport properties of materials from first principles. It solves the **Peierls-Boltzmann Transport Equation (BTE)** for phonons, moving beyond the relaxation time approximation (RTA) to capture full scattering processes, including three-phonon interactions. It is particularly adept at handling multiscale systems, from bulk crystals to thin films, superlattices, and alloys.

**Scientific domain**: Phonon Transport, Thermal Conductivity, Nanoscale Heat Transfer
**Target user community**: Researchers in thermoelectrics, thermal management, and phononics

## Theoretical Methods
- **Phonon BTE**: Solves the linearized BTE for the phonon deviation distribution function.
- **Iterative Solver**: Self-consistent solution to capture Normal (N) and Umklapp (U) scattering processes.
- **Second & Third-Order Force Constants**: Inputs derived from DFT (via codes like VASP, Quantum ESPRESSO, or Phono3py).
- **Kappa Decomposition**: Resolves thermal conductivity by mode, mean free path, and frequency.

## Capabilities
- **Bulk Properties**:
  - Lattice thermal conductivity tensor ($\kappa_{\alpha\beta}$).
  - Cumulative $\kappa$ vs. MFP.
  - Phonon lifetimes and line widths.
- **Nanostructures**:
  - Thin films (Cross-plane and In-plane).
  - Superlattices (effective medium or explicit).
  - Nanowires (diffuse boundary scattering).
- **Disorder**:
  - Virtual Crystal Approximation (VCA) for mass disorder (isotopes/alloys).

## Key Strengths
- **Beyond RTA**: One of the few publicly available codes that rigorously solves the full scattering matrix inversions.
- **Multiscale**: Can treat ballistic-diffusive regimes in nanostructures.
- **Interface**: User-friendly Python/XML inputs and HDF5 outputs.
- **Efficiency**: Parallelized with MPI and OpenMP for large q-grids.

## Inputs & Outputs
- **Inputs**:
  - `FORCE_CONSTANTS_2ND` and `FORCE_CONSTANTS_3RD` (from Phono3py/ShengBTE format).
  - XML input file defining the grid and temperature.
- **Outputs**:
  - HDF5 files containing mode-resolved properties.
  - Text summaries of thermal conductivity.

## Interfaces & Ecosystem
- **Upstream**:
  - **Phono3py**: Common generator for the necessary force constants.
  - **VASP/QE**: Source of force calculations.
- **Downstream**:
  - Plotting scripts included for visualization.

## Performance Characteristics
- **Computational Cost**: Dominant cost is the input DFPT/Force calculations. The BTE solution itself is relatively fast (minutes to hours).
- **Scaling**: Good MPI scaling for the BTE solver.

## Comparison with Other Codes
- **vs. ShengBTE**: almaBTE is considered a successor or alternative, with better object-oriented design (C++) and more features for nanostructures (thin films).
- **vs. Phono3py**: Phono3py focuses on the force constants and simple BTE; almaBTE offers more flexible BTE solvers (e.g., for superlattices) on top of those constants.

## Application Areas
- **Thermoelectrics**: Designing materials with low $\kappa_L$ to high $zT$.
- **Thermal Management**: Understanding heat flow in semiconductor thin films.
- **Isotope Engineering**: Effect of isotopic purity on diamond/Si thermal conductivity.

## Community and Support
- **Development**: Developed by Jesus Carrete, Bjorn Vermeersch, and collaborators.
- **Source**: Bitbucket.

## Verification & Sources
- **Repository**: [https://bitbucket.org/sousaw/almabte](https://bitbucket.org/sousaw/almabte)
- **Primary Publication**: Carrete et al., Comp. Phys. Comm. 220, 351 (2017).
- **Verification status**: ✅ VERIFIED
  - Active and widely cited.
  - Considered a standard tool in the field.
