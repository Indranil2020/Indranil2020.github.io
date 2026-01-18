## Official Resources
- **Repository**: https://github.com/K4ys4r/WannierPy
- **Related Project (wannierpy)**: https://github.com/henriquemiranda/wannierpy
- **License**: MIT License
- **Developers**: Community-driven (various variants exist).

## Overview
**WannierPy** represents a collection of Python-based tools and scripts designed to facilitate the post-processing and analysis of **Wannier90** output. While multiple forks and versions exist (e.g., `K4ys4r/WannierPy`, `henriquemiranda/wannierpy`), they share a common goal: providing a lightweight, scriptable interface to handle Hamiltonian matrices, plot band structures, and compute transport properties without the overhead of heavy compiled codes.

**Scientific domain**: Condensed matter physics, tight-binding analysis.
**Target user community**: Students and researchers needing quick, custom Python analysis of Wannier90 data.

## Theoretical Methods
- **Tight-Binding Interpolation**:
  - Fourier transform of real-space Hamiltonian $H(R)$ to reciprocal space $H(k)$.
  - Band structure diagonalization.
- **Transport**:
  - Calculation of basic conductivities (Boltzmann transport in relaxation time approximation).
  - Carrier lifetime interpolation.
- **Berry Phase**:
  - Evaluation of Berry curvature and anomalous Hall conductivity (in some variants).

## Capabilities
- **I/O Operations**: Reads standard `_hr.dat`, `_centres.xyz`, and `.eig` files.
- **Band Structure**:
  - Fast interpolation along high-symmetry paths.
  - Orbital-projected bands (fatbands).
- **Visualization**:
  - Matplotlib-based plotting of bands and DOS.
  - Visualization of Fermi surfaces.
- **Topological Analysis**: (Available in advanced forks) Calculation of Chern numbers and Berry curvature.

## Key Strengths
- **Lightweight**: Pure Python/NumPy implementation; easy to install and modify.
- **Educational**: excellent resource for learning how Wannier interpolation works "under the hood".
- **Flexibility**: Users can easily add custom observables or analysis routines directly in Python.

## Inputs & Outputs
- **Inputs**:
  - `wannier90_hr.dat` (Real-space Hamiltonian).
  - `wannier90.win` (Input file for k-path parsing).
- **Outputs**:
  - Matplotlib figures (PDF/PNG).
  - Text files with interpolated eigenvalues.

## Interfaces & Ecosystem
- **Wannier90**: The primary source of input data.
- **NumPy/SciPy**: Heavy reliance for matrix operations.
- **Matplotlib**: Used for all plotting functionalities.

## Computational Cost
- **Low**: Extremely efficient for standard unit cells; cost scales with the number of Wannier functions ($N_W^3$) and k-points.

## Comparison with Other Codes
- **vs [WannierBerri](file:///home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/TightBinding/4.1_Wannier_Ecosystem/WannierBerri.md)**: WannierBerri is a highly optimized, feature-rich production code for high-performance transport/optical calculations; WannierPy is a simpler, lightweight toolkit for basic analysis and quick plotting.
- **vs [PyWannier90](file:///home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/TightBinding/4.1_Wannier_Ecosystem/PyWannier90.md)**: PyWannier90 often refers to the official interface; WannierPy usually denotes community scripts for post-processing.

## Application Areas
- **Quick Plotting**: Rapidly checking band structures after a Wannier90 run.
- **Teaching**: Demonstrating TB interpolation concepts.
- **Custom Analysis**: Prototyping new observables before implementing them in Fortran/C++.

## Verification & Sources
- **Primary Source**: [GitHub Repository (K4ys4r)](https://github.com/K4ys4r/WannierPy)
- **Secondary Source**: [GitHub Repository (henriquemiranda)](https://github.com/henriquemiranda/wannierpy)
- **Verification Status**: ✅ VERIFIED (Community code).
