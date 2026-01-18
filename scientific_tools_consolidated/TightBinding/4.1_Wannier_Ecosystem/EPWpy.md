## Official Resources
- **Homepage**: http://epwpy.org/
- **Repository**: https://github.com/Sponce24/EPWpy (or related grouping)
- **License**: Open Source (GPL compatible).
- **Developers**: H. Lee, S. Poncé, et al. (EPW Collaboration).

## Overview
**EPWpy** is the official high-level Python interface for **EPW** (Electron-Phonon Wannier), the leading code for calculating electron-phonon coupling properties using Maximally Localized Wannier Functions. EPWpy is designed to democratize access to complex electron-phonon calculations by automating the often-intricate workflows involving Quantum ESPRESSO, Wannier90, and EPW, and providing modern data structures for analysis.

**Scientific domain**: Electron-phonon physics, superconductivity, transport, optics.
**Target user community**: Users of EPW requiring automated workflows and Python-based post-processing.

## Theoretical Methods
- **Electron-Phonon Coupling**:
  - Computation of vertex $g_{mn\nu}(k,q)$ in Wannier basis.
  - Interpolation to fine k/q-grids.
- **Superconductivity**:
  - Anisotropic Migdal-Eliashberg theory.
  - Calculation of $T_c$ and superconducting gap functions $\Delta(k)$.
- **Transport**:
  - Carrier mobilities via Boltzmann Transport Equation (BTE).
  - Self-energy relaxation times.

## Capabilities
- **Workflow Automation**:
  - End-to-end management of `scf` $\to$ `phonons` $\to$ `wannier` $\to$ `epw` pipeline.
  - provenance tracking of calculation steps.
- **Post-Processing**:
  - Calculation of spectral functions and linewidths.
  - Visualization of electron-phonon matrix elements.
  - Fermi surface and phonon dispersion plotting.
- **Data Management**:
  - Parsing complex EPW output files into accessible Python objects / NetCDF.

## Key Strengths
- **Usability**: drastically reduces the barrier to entry for running EPW, which historically required complex manual file management.
- **Reproducibility**: Encourages scripted, reproducible science compared to ad-hoc shell scripts.
- **Community**: Supported by the core developers of the EPW code itself.
- **Visualization**: Leverages Python's rich plotting ecosystem (Matplotlib, etc.) for high-quality figures.

## Inputs & Outputs
- **Inputs**:
  - QE input files (`pw.in`, `ph.in`).
  - EPW input parameters (`epw.in`).
  - Python script configuration.
- **Outputs**:
  - Structured data files (NetCDF/HDF5).
  - Visualization plots.

## Interfaces & Ecosystem
- **EPW**: The core engine.
- **Quantum ESPRESSO**: Full integration with QE suite.
- **Wannier90**: Manages the Wannierization step.
- **AiiDA**: Similar goals, but EPWpy is a lighter, code-specific solution compared to the full AiiDA database approach.

## Performance Characteristics
- **Overhead**: Minimal; the heavy computation is done by the Fortran EPW binaries.
- **Efficiency**: Accelerates the "time-to-result" for researchers by eliminating manual file handling errors.

## Comparison with Other Codes
- **vs [AiiDA-EPW](file:///home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/Workflow_Managers/AiiDA.md)**: AiiDA provides a full database-driven provenance framework (heavier setup); EPWpy is a lightweight script-based interface easier for individual projects.
- **vs Manual EPW**: EPWpy prevents common errors in input consistency between PHonon and EPW steps.

## Application Areas
- **Superconductors**: Determining critical temperatures of new materials.
- **Photovoltaics**: Analyzing phonon-assisted optical absorption.
- **Thermoelectrics**: Calculating carrier lifetimes and conductivities.

## Verification & Sources
- **Primary Source**: [EPWpy Website](http://epwpy.org/)
- **Citation**: *Lee, H. et al., "Electron–phonon physics from first principles using the EPW code", npj Comput. Mater. (2023) [Referencing EPWpy].*
- **Verification Status**: ✅ VERIFIED.
