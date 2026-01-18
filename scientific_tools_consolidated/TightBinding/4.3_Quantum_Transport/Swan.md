# Swan

## Official Resources
- **Repository**: https://github.com/wushidonguc/swan
- **License**: MIT License

## Overview
**Swan** is an open-source C++ software package designed for **nanoscale quantum electron transport** simulations. It employs the Non-Equilibrium Green's Function (NEGF) formalism coupled with a self-consistent Poisson solver to calculate the electrical characteristics of realistic nanodevices. Swan distinguishes itself by using **Wannier functions** as the basis set, which allows for accurate material representation while maintaining computational efficiency for large-scale atomistic simulations.

**Scientific domain**: Nanoelectronics, Device Physics, Quantum Transport
**Target user community**: Device engineers and physicists modeling transistors and quantum structures

## Theoretical Methods
- **NEGF Formalism**: Solves the Keldysh Green's functions for open quantum systems.
- **Poisson Equation**: Self-consistently solves for the electrostatic potential using a finite difference or finite element scheme.
- **Wannier Basis**: Uses Maximally Localized Wannier Functions (MLWFs) from Wannier90 to construct the device Hamiltonian.
- **Schrödinger-Poisson**: Iterative solution loop until convergence of charge density and potential.

## Capabilities
- **Device Simulation**:
  - I-V characteristics of FETs (FinFET, Nanowire FET, TFET).
  - Charge density profiles.
  - Band diagrams under bias.
- **Geometries**:
  - 1D Nanowires.
  - 2D Ultrathin bodies / Ribbons.
- **Materials**:
  - Silicon, Germanium, III-V semiconductors.
  - Transition metal dichalcogenides (TMDs).

## Key Strengths
- **Efficiency**: The use of Wannier functions provides a minimal basis set compared to plane waves or large Gaussian bases, enabling the simulation of larger devices (thousands of atoms).
- **Parallelism**: Efficient MPI parallelization for energy integration and bias points.
- **Modularity**: Object-oriented C++ design facilitates extension.

## Inputs & Outputs
- **Inputs**:
  - Hamiltonian Files (`_hr.dat`, `_xyz.dat` from Wannier90).
  - Device configuration file (JSON/Input script).
- **Outputs**:
  - `.dat` files for current, density, and potential.
  - Visualization files for ParaView (VTK).

## Interfaces & Ecosystem
- **Upstream**:
  - **Wannier90**: Essential for generating the material parameters.
  - **DFT Codes**: VASP, QE, etc., via Wannier90.
- **Downstream**:
  - **ParaView**: For 3D visualization of scalar fields (potential, density).

## Performance Characteristics
- **Speed**: Optimized dense/sparse matrix operations.
- **Scaling**: Scales up to hundreds of cores for energy points.

## Comparisons with Other Codes
- **vs. NanoTCAD ViDES**: ViDES is more comprehensive (includes drift-diffusion, Python interface); Swan is a focused C++ NEGF solver.
- **vs. OMEN**: OMEN is a high-performance HPC code; Swan is lighter and easier to deploy for smaller clusters.

## Community and Support
- **Development**: Maintained by researchers (check GitHub contributors).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/wushidonguc/swan](https://github.com/wushidonguc/swan)
- **Verification status**: ✅ VERIFIED
  - Active codebase.
  - Implements standard NEGF-Poisson loop.
