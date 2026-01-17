# Opticx

## Official Resources
- **Source Code**: https://github.com/xatu-code/opticx
- **Documentation**: [manual.pdf](https://github.com/xatu-code/opticx/blob/main/manual.pdf) (in repo)
- **License**: See Repository

## Overview
**Opticx** is a standalone Fortran code designed to evaluate **n-order optical conductivities** in periodic materials modeled with a **Tight-Binding** description. It serves as a companion tool to the Xatu code, enabling the calculation of non-linear optical properties and interfacing with Xatu to include **excitonic effects** in the optical response.

## Theoretical Methods
- **Optical Conductivity**: Calculation of linear and non-linear (n-order) optical response.
- **Tight-Binding**: Operates on tight-binding Hamiltonian models.
- **Wannier Interpolation**: Uses wannierized band structures for high-resolution integration.
- **Excitonic Effects**: Interfaces with Xatu to include scaling/shifting based on BSE results.

## Capabilities
- **Band Structure Input**: Reads wannierized Hamiltonians (e.g., from Wannier90).
- **Non-Linear Optics**: Calculation of high-order optical conductivities (SHG, THG, etc.).
- **Exciton Integration**: Can use output from Xatu to refine optical spectra with excitonic corrections.
- **Periodic Systems**: Optimized for crystals and 2D materials.

## Implementation & Tech Stack
- **Language**: Fortran.
- **Structure**: Standalone binary with input/output file interfaces.
- **Integration**: Designed to work in a workflow with **Xatu** and **Wannier90**.

## Citation
> J.J. Esteve-Paredes, M. A. García-Blázquez, A. J. Uría-Álvarez, M. Camarasa-Gómez and J. J. Palacios, npj Computational Materials 11, 13 (2025).
