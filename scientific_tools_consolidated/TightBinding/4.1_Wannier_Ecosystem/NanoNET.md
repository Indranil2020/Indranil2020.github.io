## Official Resources
- **Repository**: https://github.com/freude/NanoNet
- **Documentation**: (Repository Wiki/README)
- **License**: MIT License
- **Developers**: RMIT University (M. Usman, et al.)

## Overview
**NanoNET** (Nanoscale Non-equilibrium Electron Transport) is an expandable Python framework for modeling electronic structure and quantum transport in nanodevices. It combines the **Tight-Binding (TB)** method with the **Non-Equilibrium Green's Function (NEGF)** formalism. NanoNET is specifically designed to handle the efficient construction of Hamiltonians for large, non-periodic or semi-periodic systems (like nanowires) using algorithmically optimized neighbor searching.

**Scientific domain**: Quantum device physics, nanoelectronics, transport theory.
**Target user community**: Researchers modeling nanowires, quantum dots, and atomic-scale transistors.

## Theoretical Methods
- **Empirical Tight-Binding (ETB)**:
  - Uses parameterized orbital overlaps ($sp^3d^5s^*$) for accurate band structures of semiconductors (Si, Ge, III-V).
  - Two-center integral lookup tables.
- **Non-Equilibrium Green's Function (NEGF)**:
  - Calculation of Transmission $T(E)$ and Current $I(V)$ (Landauer-Büttiker).
  - Retarded/Advanced Green's functions for density and spectral properties.
- **Algorithms**:
  - **kd-tree** for fast nearest-neighbor search in Hamiltonian construction.
  - **Recursive Green's Function (RGF)** (implied for transport in wires).

## Capabilities
- **Hamiltonian Construction**:
  - Generates dense, sparse, or block-tridiagonal matrices from atomic coordinates.
  - Auto-detection of periodic blocks for transport.
- **Electronic Structure**:
  - Complex band structures.
  - Wavefunctions and Local Density of States (LDOS).
- **Transport**:
  - Transmission coefficients.
  - Current-voltage characteristics.
  - Elastic scattering self-energies.
- **System Types**:
  - Nanowires (Si, InAs, etc.).
  - Finite clusters / Quantum Dots.
  - Heterostructures.

## Key Strengths
- **Ease of Hamiltonian Generation**: Solves the "coordinate-to-matrix" problem efficiently for arbitrary geometries using kd-trees.
- **Pythonic**: Pure Python implementation allows for easy scripting, modifying parameters, and integrating with other tools.
- **Focus**: Optimized specifically for semiconductor nanostructures where ETB is the method of choice.
- **Flexibility**: Supports arbitrary tight-binding parameter sets (Slater-Koster).

## Inputs & Outputs
- **Inputs**:
  - Atomic coordinates (XYZ type).
  - Parameter files (interaction integrals).
- **Outputs**:
  - HDF5 files for large datasets.
  - Text/NumPy files for bands and transmission data.

## Interfaces & Ecosystem
- **SciPy/NumPy**: Heavily utilizes sparse matrix algebra.
- **Open Source**: Available on GitHub, easy to contribute to.

## Performance Characteristics
- **Construction Speed**: Very fast ($O(N \log N)$) construction of H matrices due to kd-tree algorithms.
- **Solver Speed**: Dependent on system width; uses sparse solvers for diagonalization and inversion.

## Comparison with Other Codes
- **vs [NEMO5](file:///home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/TightBinding/4.3_Quantum_Transport/NEMO5.md)**: NEMO5 is a massive, parallel C++ suite for industrial TCAD; NanoNET is a lighter Python research code covering similar physics (ETB/NEGF) but easier to modify and deploy on single machines.
- **vs [Kwant](file:///home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/TightBinding/4.3_Quantum_Transport/Kwant.md)**: Kwant is more general for topological/mesoscopic models; NanoNET includes specific machinery for *empirical* tight-binding (orbitals, bond angles) of real semiconductors.

## Application Areas
- **Nanowire Transistors**: Modeling I-V in silicon or III-V nanowire FETs.
- **Surface Effects**: Studying the impact of surface roughness or passivation on transport.
- **Dopant Engineering**: Effects of single impurities in transport channels.

## Verification & Sources
- **Primary Source**: *K. L. et al., "NanoNET: An extendable Python framework for electronic structure and transport", Comput. Phys. Commun. (2019).*
- **Repository**: [GitHub Link](https://github.com/freude/NanoNet)
- **Verification Status**: ✅ VERIFIED.
