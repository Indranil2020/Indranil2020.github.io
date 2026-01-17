# PhotoionizationGTO.jl

## Official Resources
- **Repository**: https://github.com/antoine-levitt/PhotoionizationGTO.jl
- **License**: MIT (implied by Julia ecosystem norms, check repo)
- **Documentation**: GitHub README

## Overview
PhotoionizationGTO.jl is a specialized Julia package for calculating photoionization spectra using Time-Dependent Density Functional Theory (TDDFT). It uniquely employs Gaussian-type orbitals (GTOs) to describe the bound states while handling the continuum states appropriate for photoionization processes, interfacing with the PySCF library for integrals.

**Scientific domain**: Molecular photoionization, chemical physics
**Target user community**: Julia developers, quantum chemists studying photo-ejection

## Theoretical Methods
- **Linear-Response TDDFT**: Frequency-domain response formalism
- **Gaussian Basis Sets**: Use of GTOs for efficient molecular description
- **Stieltjes Imaging / Complex Basis**: Techniques for extracting continuum information from L^2 basis sets (check specific implementation details in source)
- **Dyson Orbitals**: Calculation of Dyson orbitals for ionization amplitudes

## Capabilities
- **Photoionization Cross-Sections**: Calculation of energy-dependent cross-sections
- **Dyson Orbitals**: Analysis of the ionized electron channel
- **Julia-Python Hybrid**: Leverages Julia's speed and PySCF's robust integrals
- **Reproducibility**: Explicit dependency management via Manifest.toml

## Inputs & Outputs
- **Input formats**: Julia scripts (`.jl`), integrating PySCF molecule definitions
- **Output data types**:
  - Cross-section data arrays
  - Dyson orbital coefficients
  - Spectral plots (via Julia plotting libraries)

## Interfaces & Ecosystem
- **PySCF**: Critical dependency for ground state SCF and integrals (via PyCall)
- **Julia Ecosystem**: Integration with standard Julia linear algebra and plotting tools

## Performance Characteristics
- **Efficiency**: Julia provides near-C performance for the custom TDDFT logic.
- **Bottlenecks**: Integral evaluation depends on PySCF (C/Python) performance.

## Usage & Best Practices
- **Installation**: Use Julia's Pkg mode: `] dev https://github.com/antoine-levitt/PhotoionizationGTO.jl`
- **Workflow**: Define molecule in PySCF, pass to PhotoionizationGTO, run spectral calculation.

## Limitations & Known Constraints
- **Documentation**: Limited to README; requires reading source code for advanced features.
- **Maturity**: Research code, may lack convenience features of established packages.

## Citations
- **Primary**: Cite the repository and A. Levitt et al. publications related to the method.
