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

## Advanced Features

### Space-Dependent BTE:
- Finite Volume Method (FVM) solver
- Monte Carlo particle tracking
- Anisotropic Mean-Free-Path (aMFP) formulation
- Temperature and heat flux field calculations

### Geometry Handling:
- Complex 3D geometries via mesh files
- Porous media and phononic crystals
- Boundary condition specification
- Device-level thermal modeling

### Physics Modeling:
- Ballistic-diffusive crossover
- Boundary scattering (diffuse/specular)
- Size effects on thermal conductivity
- Interface thermal resistance

### GPU Acceleration:
- PyTorch-based linear solvers
- Efficient for large mesh calculations
- Scalable to complex geometries

## Performance Characteristics
- **Speed**: Efficient aMFP formulation makes it feasible for 3D meshes.
- **Scaling**: Scales with mesh size ($N_{vol}$) and number of phonon modes.
- **GPU support**: PyTorch acceleration available
- **Memory**: Depends on mesh resolution

## Computational Cost
- **Mesh generation**: Preprocessing step
- **BTE solution**: Minutes to hours for 3D
- **Visualization**: Fast with VTK/ParaView
- **Overall**: Efficient for device-level simulations

## Comparison with Other Codes
- **vs. almaBTE**: Both solve space-dependent BTE; OpenBTE emphasizes the FVM/aMFP approach and Python integration, while almaBTE uses Monte Carlo.
- **vs. ShengBTE**: ShengBTE is for bulk material properties; OpenBTE takes those properties and applies them to specific device geometries.
- **Unique strength**: Device-level thermal modeling with complex geometries

## Best Practices

### Workflow:
- Start with bulk BTE calculation (ShengBTE/phono3py)
- Generate appropriate mesh for geometry
- Define boundary conditions carefully
- Validate with analytical solutions when possible

### Mesh Design:
- Use appropriate mesh resolution
- Refine mesh in critical regions
- Check mesh convergence
- Balance accuracy vs computational cost

## Application Areas
- Device-level thermal modeling
- Nanostructured thermoelectrics
- Phononic crystal design
- Thermal interface engineering
- Heat spreader optimization

## Community and Support
- **Development**: MIT / UIUC (Giuseppe Romano)
- **License**: GPL-3.0
- **Repository**: GitHub (active)
- **Documentation**: https://openbte.readthedocs.io/
- **Support**: GitHub issues
- **User base**: Nanoscale thermal transport community
- **Integration**: Python ecosystem

## Verification & Sources
- **Repository**: [https://github.com/OpenBTE/OpenBTE](https://github.com/OpenBTE/OpenBTE)
- **Primary Publication**: G. Romano et al., Phys. Rev. B (2021).
- **Verification status**: ✅ VERIFIED
  - Active development and documentation.
