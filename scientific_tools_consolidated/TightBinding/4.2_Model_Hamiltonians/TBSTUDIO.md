# TBSTUDIO (Tight-Binding Studio)

## Official Resources
- **Homepage**: https://tight-binding.com/
- **Repository**: https://github.com/mohammadnakhaee/tbstudio
- **License**: Academic/Commercial (See Website)

## Overview
**TBSTUDIO** is a comprehensive software package centered around a graphical user interface (GUI) for the construction of tight-binding Hamiltonians from first-principles data. It simplifies the often complex workflow of fitting Slater-Koster parameters to Density Functional Theory (DFT) band structures. By automating the fitting process and providing visualization tools, it serves as a bridge between ab initio codes (like VASP, Quantum ESPRESSO) and model analysis tools.

**Scientific domain**: Materials Modeling, Band Structure Fitting
**Target user community**: Materials scientists needing to extract effective models from DFT

## Theoretical Methods
- **Slater-Koster Formalism**: Parameterization of hopping integrals based on orbital symmetries ($ss\sigma$, $pp\pi$, etc.).
- **Fitting Algorithm**: Levenberg-Marquardt non-linear least squares optimization to minimize the difference between DFT and TB eigenvalues.
- **Basis Sets**: Support for orthogonal and non-orthogonal bases.
- **Spin-Orbit Coupling**: Inclusion of atomic SOC parameters.

## Capabilities
- **Model Generation**:
  - Auto-fit bands from VASP/QE/Wien2k.
  - Generates Slater-Koster hopping tables.
- **Visualization**:
  - 3D rendering of the crystal structure and orbital positions.
  - Interactive plots comparing DFT and fitted TB bands.
- **Exports**:
  - Python (Pybinding compatible output).
  - MATLAB, C++, Fortran, Mathematica.
  - Raw Hamiltonian matrices.

## Key Strengths
- **GUI-Driven**: Makes the sophisticated task of tight-binding parameterization accessible to users without deep coding experience.
- **Cross-Platform**: Runs on Windows, Linux, and macOS.
- **Integration**: Designed to feed into widely used solvers like Pybinding or Green's function codes.

## Inputs & Outputs
- **Inputs**:
  - Crystal structure (POSCAR, CIF).
  - Band structure data (EIGENVAL).
- **Outputs**:
  - Fitted parameter files.
  - Source code defining the model in various languages.

## Interfaces & Ecosystem
- **Upstream**: VASP, Quantum ESPRESSO, Abinit, Wien2k.
- **Downstream**: Pybinding, Kwant (via script), custom codes.

## Performance Characteristics
- **Efficiency**: Fitting is performed locally; speed depends on the number of orbitals and k-points fitted.
- **Usability**: Interactive feedback loop significantly speeds up the model generation process compared to command-line fitting tools.

## Comparison with Other Codes
- **vs. Wannier90**: Wannier90 is the gold standard for exact (interpolated) tight-binding models. TBSTUDIO uses the Slater-Koster *approximation*, which is less exact but physically more intuitive (fewer, longer-range parameters) and often transfers better to varying geometries.
- **vs. PythTB**: PythTB is a library for *using* models. TBSTUDIO is a tool for *creating* them.

## Application Areas
- **Heterostructures**: Creating transferrable models for interface calculations.
- **Device Simulation**: Generating input Hamiltonians for transport codes.

## Community and Support
- **Development**: Mohammad Nakhaee.
- **Source**: GitHub / Website.

## Verification & Sources
- **Website**: [https://tight-binding.com/](https://tight-binding.com/)
- **Primary Publication**: M. Nakhaee et al., arXiv:1910.02917.
- **Verification status**: ✅ VERIFIED
  - Active tool with commercial/academic licensing dual model.
