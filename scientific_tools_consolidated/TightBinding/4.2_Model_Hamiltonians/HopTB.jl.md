# HopTB.jl

## Official Resources
- **Homepage**: https://hoptb.github.io/HopTB.jl/dev/
- **Repository**: https://github.com/HopTB/HopTB.jl
- **License**: MIT License

## Overview
**HopTB.jl** is a Julia package designed for constructing and analyzing tight-binding Hamiltonians, with a unique focus on **non-orthogonal bases**. It serves as a bridge between first-principles Density Functional Theory (DFT) codes and model physics, allowing users to import Hamiltonians from **Wannier90**, **OpenMX**, and **FHI-aims**. Beyond standard band structures, HopTB.jl provides a powerful suite of tools for calculating linear and **non-linear response functions**, including optical conductivity, Hall effects, and second harmonic generation.

**Scientific domain**: Theoretical Materials Science, Non-Linear Optics
**Target user community**: Researchers bridging DFT and effective model calculations

## Theoretical Methods
- **Tight-Binding**: Supports both standard orthogonal ($S_{ij} = \delta_{ij}$) and generalized non-orthogonal ($S_{ij} \neq 0$) tight-binding models.
- **Linear Response**: Kubo formula implementation for conductivity tensors.
- **Geometric Phase**: Calculation of Berry curvature, Berry connection, and their dipole moments.
- **Non-Linear Optics**:
  - Shift Current (bulk photovoltaic effect).
  - Second Harmonic Generation (SHG).
  - Injection Current.

## Capabilities
- **Interfaces**: 
  - Read generic Wannier90 outputs (`_hr.dat`, `.win`).
  - Read OpenMX and FHI-aims tight-binding formats.
- **Observables**:
  - Band structures, Fermi surfaces.
  - Optical Conductivity $\sigma(\omega)$.
  - Anomalous Hall Continuity.
  - Spin Hall Conductivity.
  - Non-linear conductivities (Shift, Berry Dipole).
- **Symmetrization**: Tools to enforce crystal symmetries on tight-binding models derived from numerical data.

## Key Strengths
- **Non-Orthogonality**: One of the few transport/response codes that correctly handles the overlap matrix $S$ from local-orbital DFT codes (OpenMX/FHI-aims) without assuming orthogonality.
- **Non-Linear Optics**: Specialized features for the emerging field of non-linear Hall effects and shift currents, which are not found in standard TB packages.
- **Julia Performance**: Exploits Julia's JIT compilation for efficient integration over dense k-grids.

## Inputs & Outputs
- **Inputs**:
  - DFT output files (Wannier/OpenMX/FHI-aims).
  - Julia scripts defining calculation parameters.
- **Outputs**:
  - Computed tensors (conductivity, etc.).
  - Plotting objects.

## Interfaces & Ecosystem
- **Upstream**: Wannier90, OpenMX, FHI-aims.
- **Downstream**: Julia plotting libraries (Plots.jl, Makie.jl).

## Performance Characteristics
- **Efficiency**: Highly optimized for k-point summation.
- **Parallelism**: Julia multi-threading support.

## Comparison with Other Codes
- **vs. WannierBerri**: WannierBerri is the gold standard for Berry phase properties from Wannier90. HopTB.jl offers similar capabilities around optics and Hall effects but adds support for non-orthogonal bases (non-Wannier90 sources).
- **vs. TB2J**: TB2J calculates magnetic parameters; HopTB.jl calculates optical/electronic response.

## Application Areas
- **Topological Photovoltaics**: Studying shift currents in Weyl semimetals.
- **Spintronics**: Spin Hall effect in complex oxides.
- **Methodology**: Testing the validity of orthogonal approximations in tight-binding.

## Community and Support
- **Development**: HopTB team (GitHub).
- **Source**: GitHub.

## Verification & Sources
- **Documentation**: [https://hoptb.github.io/HopTB.jl/dev/](https://hoptb.github.io/HopTB.jl/dev/)
- **Verification status**: ✅ VERIFIED
  - Active Julia package.
