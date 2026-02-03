# OpenBTE

## Official Resources
- **Homepage**: https://github.com/OpenBTE/OpenBTE
- **Documentation**: https://openbte.readthedocs.io/
- **License**: GPL-3.0

## Overview
**OpenBTE** is an open-source vibrational transport solver designed to compute **lattice thermal conductivity** and heat transport maps in **multidimensional nanostructures**. Unlike bulk BTE solvers (like ShengBTE), OpenBTE solves the **space-dependent Boltzmann Transport Equation** for phonons, making it capable of modeling size effects, boundary scattering, and heat flow tailored geometries (membranes, porous materials, nanowires).

**Scientific domain**: Nanoscale Heat Transport, Thermoelectrics
**Target user community**: Researchers bridging material properties and device geometry

## Theoretical Methods
- **Space-Dependent BTE**: Solves $\mathbf{v} \cdot \nabla T + \frac{T - T_0}{\tau} = 0$ (in RTA) or more complex forms.
- **Solvers**:
  - **Finite Volume Method (FVM)**: For deteminstic solution on a mesh.
  - **Monte Carlo**: For particle-based tracking.
- **aMFP**: Anisotropic Mean-Free-Path formulation to reduce angular variables.

## Capabilities
- **Observables**:
  - Effective Thermal Conductivity ($\kappa_{eff}$).
  - Temperature maps $T(\mathbf{r})$.
  - Heat flux fields $\mathbf{J}(\mathbf{r})$.
- **Geometries**:
  - 1D/2D/3D complex shapes (defined by meshes).
  - Porous media (phononic crystals).
- **Physics**:
  - Boundary scattering (diffuse/specular).
  - Ballistic-to-diffusive crossover.

## Key Strengths
- **Geometry Awareness**: Can simulate real device shapes, not just bulk unit cells.
- **Ab Initio Link**: Directly uses phonon lifetimes/velocities from `ShengBTE` or `Phono3py` as material inputs.
- **Optimization**: GPU acceleration via PyTorch for linear solvers.

## Inputs & Outputs
- **Inputs**:
  - Bulk phonon properties (BTE solution for bulk).
  - Mesh files (`.msh`).
- **Outputs**:
  - HDF5/VTK files for visualization in ParaView.

## Interfaces & Ecosystem
- **Upstream**: ShengBTE, Phono3py.
- **Python**: Fully Pythonic API.

## Performance Characteristics
- **Speed**: Efficient aMFP formulation makes it feasible for 3D meshes.
- **Scaling**: Scales with mesh size ($N_{vol}$) and number of phonon modes.

## Comparison with Other Codes
- **vs. almaBTE**: Both solve space-dependent BTE; OpenBTE emphasizes the FVM/aMFP approach and Python integration, while almaBTE uses Monte Carlo.
- **vs. ShengBTE**: ShengBTE is for bulk material properties; OpenBTE takes those properties and applies them to specific device geometries.

## Community and Support
- **Development**: MIT / UIUC (Giuseppe Romano).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/OpenBTE/OpenBTE](https://github.com/OpenBTE/OpenBTE)
- **Primary Publication**: G. Romano et al., Phys. Rev. B (2021).
- **Verification status**: ✅ VERIFIED
  - Active development and documentation.
