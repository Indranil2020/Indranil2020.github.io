# fem-tddft

## Official Resources
- Source Repository: https://github.com/brandonwood/fem-tddft
- Documentation: https://github.com/brandonwood/fem-tddft (README)
- License: GNU General Public License v3.0

## Overview
fem-tddft is a real-space Density Functional Theory (DFT) and Time-Dependent DFT (TDDFT) code built upon the finite-element method (FEM). It leverages the deal.II finite element library to provide high spatial resolution and flexible mesh adaptivity, making it particularly suitable for studying materials from first principles where complex geometries or local field enhancements are important.

**Scientific domain**: Nanoplasmonics, materials modeling, excited electron dynamics
**Target user community**: Researchers studying optical properties and electron dynamics in complex geometries

## Theoretical Methods
- real-time Time-Dependent Density Functional Theory (RT-TDDFT)
- Kohn-Sham DFT
- Finite Element Method (FEM) discretization (via deal.II)
- Real-space formulation
- All-electron capability (implied by FEM framework flexibility)

## Capabilities
- Ground state electronic structure
- Real-time electron dynamics
- Optical response calculations
- Adaptive mesh refinement
- Parallel execution (MPI)

## Key Strengths
### Finite Element Basis
- High local accuracy with adaptive mesh refinement
- Flexible handling of complex boundary conditions and geometries
- Systematic convergence

### Real-Space Dynamics
- Direct simulation of time-dependent phenomena
- Suitable for strong fields and nonlinear response

## Inputs & Outputs
- **Input**: Parameter files (prm), mesh files
- **Output**: Density, dipole moments over time, energy, visualization files (VTK/Vtu)

## Computational Cost
- **Method**: Sparse matrix operations typical of FEM.
- **Scalability**: Parallelized using MPI, scales to computing clusters.

## Interfaces & Ecosystem
- **Dependencies**: Built on **deal.II** finite element library (C++).
- **Visualization**:
  - Outputs `.vtk` and `.vtu` files for Paraview/VisIt.
- **Input**:
  - Parameter files (`.prm`) standard to deal.II applications.

## Advanced Features
- **Nanoplasmonics**:
  - Tailored for metallic nanoparticles and local field enhancements.
- **Adaptive Mesh Refinement**:
  - Dynamically refines grid in regions of high electron density gradient.

## Community and Support
- **Status**: Research code, active development by Brandon Wood et al.
- **License**: GPL v3.0 ensures open availability.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/brandonwood/fem-tddft
2. Linked publications on repository

**Confidence**: VERIFIED
**Status**: Open Source, Research code
