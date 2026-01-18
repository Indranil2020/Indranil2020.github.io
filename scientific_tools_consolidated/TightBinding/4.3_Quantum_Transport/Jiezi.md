# Jiezi

## Official Resources
- **Repository**: https://github.com/Jiezi-negf/Jiezi
- **License**: GNU General Public License v3.0

## Overview
**Jiezi** is an open-source Python software package for the **quantum transport simulation** of nanoscale devices, specifically optimized for **Gate-All-Around Carbon Nanotube Field-Effect Transistors (GAA-CNTFETs)**. It combines the Non-Equilibrium Green's Function (NEGF) method for quantum transport with a Finite Element Method (FEM) solver for the 3D Poisson equation, enabling the rigorous self-consistent simulation of electrostatics and current flow in complex 1D channel architectures.

**Scientific domain**: Cntfet Modeling, Nanoelectronics, Device Simulation
**Target user community**: Device engineers designing carbon nanotube transistors

## Theoretical Methods
- **NEGF**: Solves for the carrier density and transmission in the ballistic limit (or with scattering).
- **Mode Space Approach**: Transforms the Hamiltonian of the nanotube into decoupled modes to speed up calculations.
- **Finite Element Method (FEM)**: Solves the 3D Poisson equation for the electrostatic potential, handling complex gate geometries and dielectrics.
- **Self-Consistency**: Anderson mixing or Newton-Raphson iteration to converge charge and potential.

## Capabilities
- **Device Characteristics**:
  - $I_d-V_g$ (Transfer) and $I_d-V_d$ (Output) curves.
  - Subthreshold Swing (SS) and DIBL analysis.
  - Quantum capacitance extraction.
- **Geometry Support**:
  - Cylindrical Gate-All-Around inputs.
  - Variable oxide thicknesses and high-k dielectrics.
- **Physics**:
  - Band-to-band tunneling (BTBT).
  - Schottky barrier effects at contracts.

## Key Strengths
- **Python-Based**: Easy to understand, modify, and integrate with optimization scripts (e.g., SciPy).
- **Specialization**: tailored specifically for cylindrical geometries common in nanotubes/nanowires, offering efficiency advantages over generic solvers.
- **Meshing**: Uses `gmsh` or similar tools for flexible FE meshing of the dielectric environment.

## Inputs & Outputs
- **Inputs**:
  - Device geometry specification (channel length, diameter, oxide).
  - Bias conditions.
  - Mesh files.
- **Outputs**:
  - Current values.
  - Potential distribution maps.
  - Energy-resolved electron density.

## Interfaces & Ecosystem
- **Dependencies**: NumPy, SciPy, Matplotlib.
- **Meshing**: Compatible with standard mesh formats.

## Performance Characteristics
- **Speed**: Mode-space approach makes it significantly faster than real-space atomistic solvers for long channels.
- **Parallelism**: Energy point integration can be parallelized.

## Comparison with Other Codes
- **vs. NanoTCAD ViDES**: ViDES is more general (Graphene, MoS2, Si); Jiezi is highly focused on CNTs.
- **vs. Feldman's code**: Similar functionality but Jiezi offers a modern Python implementation.

## Application Areas
- **Digital Logic**: Assessment of CNTFETs for sub-5nm technology nodes.
- **Sensors**: Modeling the electrostatic gating effect of biomolecules on CNTs.

## Community and Support
- **Development**: Developed by academic researchers (Peking University group).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/Jiezi-negf/Jiezi](https://github.com/Jiezi-negf/Jiezi)
- **Verification status**: ✅ VERIFIED
  - Functional repository with documentation.
