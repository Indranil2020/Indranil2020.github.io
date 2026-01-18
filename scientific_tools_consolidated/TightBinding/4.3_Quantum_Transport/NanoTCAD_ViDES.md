# NanoTCAD ViDES

## Official Resources
- **Homepage**: http://vides.nanotcad.com/
- **Repository**: https://github.com/vides-hub/vides
- **License**: BSD License

## Overview
**NanoTCAD ViDES** (Vintage Integrated Development Environment for Simulations) is an open-source software package for the simulation of nanoscale electronic devices. It is particularly renowned for its ability to simulate **2D material-based devices** (graphene, MoS2) using the **Non-Equilibrium Green's Function (NEGF)** formalism self-consistently coupled with a 2D/3D **Poisson solver**. The code is wrapped in Python, providing a flexible scripting environment for investigating novel transistor architectures.

**Scientific domain**: Nanoelectronics, CNTs, Graphene, TMDs
**Target user community**: Device engineers and physicists developing post-silicon logic

## Theoretical Methods
- **NEGF Formalism**: Standard coherent transport formalism for calculating density and current.
- **Poisson Solver**: Finite difference solver for the electrostatic potential in 2D/3D geometries.
- **Self-Consistency**: Newton-Raphson or Gummel iteration to solve the coupled Schrödinger-Poisson system.
- **Hamiltonians**:
  - Tight-Binding models (Nearest neighbor, etc.).
  - Massless Dirac Fermions (Continuous models).
  - Maximally Localized Wannier Functions (integration).

## Capabilities
- **Simulations**:
  - Graphene Nanoribbon FETs (GNR-FETs).
  - Carbon Nanotube FETs (CNT-FETs).
  - TMD Transistors (MoS2, WSe2).
  - Heterojunctions and tunneling barriers.
- **Observables**:
  - Transfer characteristics ($I_d-V_g$).
  - Output characteristics ($I_d-V_d$).
  - Local Density of States (LDOS).
  - Potential profiles and subband structures.

## Key Strengths
- **Python Interface**: The `pyViDES` module allows users to construct simulations using standard Python syntax, making it highly accessible and easy to integrate with plotting libraries.
- **2D Focus**: Specialized routines and material parameters for graphene and transition metal dichalcogenides.
- **Drift-Diffusion**: Also includes a semi-classical drift-diffusion module for comparing ballistic vs classical limits.

## Inputs & Outputs
- **Inputs**: Python scripts defining the device geometry, materials, and bias loop.
- **Outputs**:
  - Text files (currents).
  - Grid data (potential, charge) for visualization.

## Interfaces & Ecosystem
- **Wannier90**: Can import Wannier Hamiltonians for atomistic accuracy.
- **Python**: Full integration with NumPy/SciPy/Matplotlib.

## Performance Characteristics
- **Computational Cost**: NEGF inversion is $O(N_y^3)$ (width). Efficient for narrow ribbons/nanotubes; slower for wide devices.
- **Parallelism**: MPI parallelization over energy points.

## Comparison with Other Codes
- **vs. Kwant**: Kwant generally lacks the built-in self-consistent Poisson solver required for realistic transistor characteristics (I-V curves); ViDES provides this "TCAD" functionality out-of-the-box.
- **vs. NEMO5**: NEMO5 is a heavier, industrial-scale code; ViDES is lighter and better suited for rapid academic prototyping of 2D devices.

## Community and Support
- **Development**: University of Pisa (Gianluca Fiori, Giuseppe Iannaccone).
- **Source**: GitHub and website.

## Verification & Sources
- **Website**: [http://vides.nanotcad.com/](http://vides.nanotcad.com/)
- **Primary Publication**: G. Fiori and G. Iannaccone, IEEE Electron Device Lett. (2007).
- **Verification status**: ✅ VERIFIED
  - Well-established in the 2D device community.
