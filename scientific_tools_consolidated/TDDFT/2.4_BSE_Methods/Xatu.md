# Xatu

## Official Resources
- **Source Code**: https://github.com/xatu-code/xatu
- **Paper**: [Efficient computation of optical excitations in two-dimensional materials with the Xatu code](https://doi.org/10.1016/j.cpc.2023.109001) (Computer Physics Communications)
- **License**: Ask authors / See Repository

## Overview
Xatu (**eXcitons from ATomistic calcUlations**) is a specialized program and library designed to solve the **Bethe-Salpeter Equation (BSE)** for solids to obtain the exciton spectrum. It is particularly optimized for **two-dimensional (2D) materials** and operates as a post-processing tool taking electronic band structures from Tight-Binding models or DFT calculations (based on local orbitals) as input.

## Theoretical Methods
- **Bethe-Salpeter Equation (BSE)**: Full solution for excitonic spectra.
- **Tight-Binding / LCAO**: Uses local orbital basis sets (from codes like SIESTA or TB models).
- **Dielectric Screening**: Implements various screening models suitable for 2D and 3D systems.
- **Optical Properties**: Calculates excitonic absorption spectra and wavefunctions.

## Capabilities
- **Exciton Spectrum**: Diagonalization of the BSE Hamiltonian.
- **Dimensionality**: Specialized for 2D materials (e.g., hBN, MoS2) but applicable to general solids.
- **Input Interfaces**: Compatible with local orbital band structures (e.g., Wannier90, Tight-Binding).
- **Analysis**: Characterization of exciton binding energies, wavefunctions, and optical strengths.

## Implementation & Tech Stack
- **Languages**: C++, Fortran, Python.
- **Libraries**: Built upon **Armadillo** C++ algebra library (BLAS, LAPACK, ARPACK backend).
- **Parallelization**: Designed for efficiency on modern architectures.

## Workflow
1.  **DFT/TB Calculation**: Generate electronic structure (Hamiltonian in local basis).
2.  **Screening**: Compute or model the static dielectric function.
3.  **BSE Construction**: Xatu builds the electron-hole interaction kernel.
4.  **Diagonalization**: Solves for exciton eigenstates and eigenvalues.
5.  **Post-processing**: Analysis of optical absorption and exciton character.

## Application Areas
- **2D Materials**: Transition metal dichalcogenides (TMDs), hBN, graphene derivatives.
- **Excitonics**: Study of bound electron-hole pairs in reduced dimensions.
- **Optical Materials**: Designing materials for optoelectronics.

## Citation
If you use Xatu, please cite:
> A. J. Uría-Álvarez et al., "Efficient computation of optical excitations in two-dimensional materials with the Xatu code", Computer Physics Communications (2023).
