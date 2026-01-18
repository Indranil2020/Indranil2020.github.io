# Spirit

## Official Resources
- **Homepage**: https://spirit-docs.readthedocs.io/
- **Repository**: https://github.com/spirit-code/spirit
- **License**: MIT License

## Overview
**Spirit** is a comprehensive, cross-platform framework for **atomistic spin dynamics** simulations. It combines a highly efficient simulation core with powerful visualization capabilities, setting a new standard for ease of use in nanomagnetism research. Spirit allows for the simulation of complex magnetic textures—such as **skyrmions**, **domain walls**, and **spin spirals**—and their dynamics under thermal fluctuations and external driving forces. Distinctively, it offers a **desktop GUI** and **web interface** for real-time interaction with the spin system.

**Scientific domain**: Nanomagnetism, Spintronics, Skyrmionics
**Target user community**: Researchers studying magnetic stability, dynamics, and phase transitions

## Theoretical Methods
- **Landau-Lifshitz-Gilbert (LLG)**: Stochastic differential equation for spin dynamics at finite temperature.
- **Geodesic Nudged Elastic Band (GNEB)**: Method for finding Minimum Energy Paths (MEP) and activation barriers between magnetic states.
- **Minimization**: Direct energy minimization (LBFGS, VP) to find ground states.
- **Hamiltonian**:
  - Heisenberg Exchange ($J_{ij}$).
  - Dzyaloshinskii-Moriya Interaction (DMI).
  - Dipole-Dipole Interaction (long-range).
  - Zeeman and Anisotropy terms.
  - Higher-order interactions (4-spin).

## Capabilities
- **Simulations**:
  - **Dynamics**: Time evolution of spins under fields, currents (STT), or thermal noise.
  - **Transitions**: Calculation of energy barriers and lifetimes of magnetic bits (skyrmions).
  - **Structure**: Finding complex ground states (frustrated magnets, spin lattices).
- **Visualization**:
  - Real-time 3D rendering of spins.
  - Interactive manipulation (poke/drag spins).
- **Systems**:
  - Thin films, multilayers, and defects.
  - Periodic boundaries or finite islands.

## Key Strengths
- **Interactivity**: The ability to interact with the simulation in real-time is unique and invaluable for building intuition about complex magnetic textures.
- **Performance**: High-performance C++ core with **GPU (CUDA)** and OpenMP parallelism.
- **Methodology**: Implements state-of-the-art GNEB and Harmonic Transition State Theory (HTST) methods for quantitative stability analysis.
- **Cross-Platform**: Runs on Linux, Windows, macOS, and even in the browser (WebAssembly).

## Inputs & Outputs
- **Inputs**:
  - Configuration files defining geometry and Hamiltonian parameters ($J$, $D$, $K$).
  - Python scripts (`import spirit`).
- **Outputs**:
  - Spin configuration files (OVF, flexible text).
  - Energy logs and reaction paths.

## Interfaces & Ecosystem
- **Python**: Full control via a Python API.
- **OVF**: Compatible with the OOMMF Vector Format standard.

## Performance Characteristics
- **Speed**: GPU acceleration enables the simulation of millions of spins with interactive frame rates.
- **Scalability**: $O(N)$ for short-range interactions; $O(N \log N)$ for dipolar interactions (FFT).

## Comparison with Other Codes
- **vs. UppASD**: UppASD focuses more on thermodynamics (phase transitions, specific heat) and material specificity; Spirit excels in energy landscapes (barriers/GNEB) and visualization.
- **vs. Mumax3**: Mumax3 is a continuum (micromagnetic) code; Spirit is atomistic (discrete spins), essential for calculating properties of very small textures like atomic-scale skyrmions.

## Community and Support
- **Development**: Forschungszentrum Jülich (Stefan Blügel group).
- **Source**: GitHub.

## Verification & Sources
- **Website**: [https://spirit-docs.readthedocs.io/](https://spirit-docs.readthedocs.io/)
- **Primary Publication**: G. P. Müller et al., Phys. Rev. B 99, 224414 (2019).
- **Verification status**: ✅ VERIFIED
  - Active, modern, and widely used tool.
