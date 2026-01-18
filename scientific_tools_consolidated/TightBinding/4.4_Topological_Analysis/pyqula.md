# pyqula

## Official Resources
- **Homepage**: https://github.com/joselado/pyqula
- **Repository**: https://github.com/joselado/pyqula
- **License**: GPL-3.0

## Overview
**pyqula** (Python Quantum Lattice) is a powerful standard-library for simulating quantum tight-binding models, with a distinctive focus on **interacting systems** and **topological phases**. Unlike many pure TB codes, pyqula includes self-consistent mean-field solvers for Hubbard and Heisenberg interactions, allowing it to simulate the interplay between correlation (magnetism, superconductivity) and topology. It also implements the **Kernel Polynomial Method (KPM)** for large-scale calculations.

**Scientific domain**: Correlated Topological Matter, Quantum Transport
**Target user community**: Theorists studying interacting topological phases (e.g., twisted bilayers, magnetic TIs)

## Theoretical Methods
- **Tight-Binding**: Slater-Koster or manual hopping.
- **Mean-Field Theory**:
  - Hubbard Model (Self-consistent $U$).
  - BdG Superconductivity (Self-consistent $\Delta$).
- **Topological Invariants**:
  - Chern Numbers (Real-space and k-space).
  - $\mathbb{Z}_2$ Invariants.
  - Local Chern Marker (for amorphous topology).
- **KPM**: Order-N expansion for DOS and Spectral Functions.

## Capabilities
- **Model Construction**:
  - Built-in generators for Honeycomb, Kagome, Square, Lieb lattices.
  - Twisted Bilayer Graphene (TBG) constructors.
- **Interactions**:
  - Non-collinear magnetism.
  - Unconventional superconductivity (d-wave, p-wave).
- **Observables**:
  - Band structures.
  - Local Density of States (LDOS).
  - Quantum Transport (Conductance).
  - Berry Curvature maps.

## Key Strengths
- **Correlations + Topology**: One of the few Python codes that seamlessly integrates mean-field order parameters into topological analysis. You can start with a TI and turn on Hubbard U to see if it becomes an antiferromagnet.
- **Real-Space Topology**: Implements the Local Chern Marker, enabling topological characterization of disordered or amorphous systems where $k$ is not a good quantum number.
- **KPM Scalability**: Can handle millions of atoms for spectral functions, overcoming the limits of exact diagonalization.

## Inputs & Outputs
- **Inputs**: Python scripts defining geometry and interaction strengths.
- **Outputs**:
  - Matplotlib plots (bands, Fermi surfaces).
  - NumPy arrays of observables.

## Interfaces & Ecosystem
- **Core**: Built on NumPy/SciPy.
- **Visualization**: Extensive internal plotting library using Matplotlib.

## Performance Characteristics
- **Speed**: Mixed. Core arithmetic is NumPy (fast), but Python loops in self-consistent cycles can be slower than Fortran. KPM is highly efficient.
- **Scalability**: KPM scales to $O(N)$; Exact Diagonalization to $O(N^3)$.

## Comparison with Other Codes
- **vs. Pybinding**: Pybinding is faster (C++) but non-interacting. Pyqula allows interactions.
- **vs. Kwant**: Kwant is better for transport in arbitrary geometries. Pyqula is better for bulk interacting phases and topological markers.

## Application Areas
- **Twisted Moiré Systems**: Simulating correlated insulating states in TBG.
- **Topological Superconductivity**: Majorana modes in hybrid nanowires.
- **Amorphous Topology**: Chern numbers in random lattices.

## Community and Support
- **Development**: Jose Lado (Aalto University).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/joselado/pyqula](https://github.com/joselado/pyqula)
- **Verification status**: ✅ VERIFIED
  - Widely cited research code.
