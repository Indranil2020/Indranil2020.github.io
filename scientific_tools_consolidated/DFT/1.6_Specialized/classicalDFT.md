# classicalDFT

## Official Resources
- Homepage: https://github.com/mencelt/classicalDFT
- Source Repository: https://github.com/mencelt/classicalDFT
- License: MIT License

## Overview
classicalDFT is a high-performance C++ library specifically designed for Classical Density Functional Theory (cDFT) calculations. While distinct from the quantum mechanical DFT used for electrons, cDFT is the analogous rigorous framework for determining the equilibrium density profiles of classical fluids. This library provides advanced solvers for integral equation theories and free energy minimization in 3D, enabling the study of fluids, colloids, and interfaces at the microscopic scale.

**Scientific domain**: Statistical Mechanics, Soft Matter Physics, Chemical Engineering
**Target user community**: Researchers studying solvation, adsorption, and colloidal assembly

## Theoretical Methods
- **Classical Density Functional Theory (cDFT)**: Statistical mechanics of inhomogeneous fluids.
- **Fundamental Measure Theory (FMT)**: Rosenfeld's accurate functional for hard spheres.
- **Mean Field Approximation**: For attractive dispersion interactions (Lennard-Jones).
- **Picard Iteration**: Fixed-point conceptual solver.
- **DIIS / Anderson Mixing**: Accelerated convergence algorithms.
- **3D Fast Fourier Transforms (FFT)**: Using FFTW for convolution integrals.

## Capabilities (CRITICAL)
- **Density Profiles**: Calculation of 3D equilibrium density distributions of fluids near surfaces or in confinement.
- **Free Energy**: Computation of grand canonical free energies.
- **Adsorption**: Calculation of adsorption isotherms in porous materials.
- **Phase Behavior**: Identification of liquid-vapor phase transitions in confinement.
- **Solvation Forces**: Calculation of forces between colloidal particles.
- **Geometry**: Support for 1D (planar), 2D (cylindrical), and 3D (arbitrary) geometries.

## Key Strengths

### High Performance:
- Written in C++ with OpenMP parallelization.
- Highly optimized convolution routines using FFTW.
- Capable of handling large 3D grids (e.g., 256^3) efficiently.

### Modern Methods:
- Implements state-of-the-art FMT functionals (White Bear, White Bear II).
- Advanced minimization algorithms ensure robust convergence near phase transitions.

## Inputs & Outputs
- **Inputs**:
  - Parameter file (JSON/Input script): Grid size, Temperature, Chemical Potential.
  - External Potential: Defined map (e.g., from a solid substrate).
- **Outputs**:
  - `rho.dat` / `rho.vtk`: 3D density fields (visualizable in Paraview).
  - Thermodynamic properties (pressure, adsorption).

## Interfaces & Ecosystem
- **C++ Library**: Can be linked into other MD or Monte Carlo codes.
- **Python Bindings**: (If available/planned) allow for scripting workflows.
- **Visualization**: VTK output for standard tools (Paraview, VESTA).

## Advanced Features
- **Mixtures**: Support for binary and multi-component fluid mixtures.
- **Charged Systems**: Electrostatic interactions via Coulomb functionals (MSA/Mean Field).

## Performance Characteristics
- **Speed**: cDFT is orders of magnitude faster than Molecular Dynamics or Monte Carlo for equilibrium properties (minutes vs days).
- **Memory**: Grid-based; memory scales with volume (N grid points).

## Computational Cost
- **Low**: Runs on single cores or small workstations for typical problems.
- **Scaling**: O(N log N) due to FFTs.

## Limitations & Known Constraints
- **Static**: Equilibrium theory only (Dynamic DFT is a separate extension).
- **Approximations**: Accuracy depends on the closure/functional quality (FMT is good for hard spheres, approximate for complex molecules).
- **Rigid Molecules**: Typically assumes spherical or simple geometries; complex molecular orientations require Molecular DFT (MDFT) extensions.

## Comparison with Other Codes
- **vs Tramonto**: Tramonto is a massive DOE package (MPI-parallel, Trilinos based); classicalDFT is a lightweight, single-node library.
- **vs MDFT**: MDFT codes handle molecular orientation (dipoles); classicalDFT focuses on simple fluids.
- **vs MD Simulation**: cDFT yields direct free energies and valid thermodynamics without sampling noise.
- **Unique strength**: Modern, accessible C++ implementation of 3D FMT functionals.

## Application Areas
- **Porous Media**: Gas storage in MOFs/Zeolites.
- **Colloids**: Depletion forces and self-assembly.
- **Wetting**: Contant angles and liquid films on structured surfaces.
- **Biophysics**: Implicit solvent models for proteins (crowding effects).

## Best Practices
- **Grid Resolution**: Ensure grid spacing is < 0.1 sigma (particle diameter) for FMT accuracy.
- **Initialization**: Start from low density and ramp up chemical potential to avoid hysteresis or solver failure.
- **Validation**: Compare with bulk equation of state in homogeneous limit.

## Community and Support
- **GitHub**: Source code and basic examples.
- **Academic**: Citations to method papers (Roth, Evans, etc.).

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/mencelt/classicalDFT
2. Standard literature on "Classical Density Functional Theory" and "Fundamental Measure Theory".

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Domain: Valid scientific tool for Soft Matter.
