# Pybinding

## Official Resources
- **Homepage**: https://docs.pybinding.site/
- **Documentation**: https://docs.pybinding.site/en/stable/
- **Repository**: https://github.com/dean0x7d/pybinding
- **License**: BSD 3-Clause License

## Overview
**Pybinding** is a high-performance Python package for numerical tight-binding calculations. It is engineered to handle **large-scale systems** (millions of atoms) by combining a user-friendly Python interface with a highly optimized C++ core. It excels at constructing arbitrary lattice geometries, applying external fields and disorder, and computing electronic properties using efficient methods like the **Kernel Polynomial Method (KPM)**.

**Scientific domain**: Mesoscopic Physics, Graphene/2D Materials
**Target user community**: Researchers simulating transport and spectra in large disordered systems

## Theoretical Methods
- **Tight-Binding Formalism**: Orthogonal tight-binding models.
- **Kernel Polynomial Method (KPM)**: Order-N expansion of spectral functions (Chebyshev polynomials) to compute Density of States (DOS) and conductivity without full diagonalization.
- **Exact Diagonalization**: Interfaces to LAPACK and ARPACK (sparse) solvers.
- **Recursion Methods**: For certain transport properties.

## Capabilities
- **Model Construction**:
  - Arbitrary crystal lattices (1D, 2D, 3D).
  - Finite and periodic systems.
  - Complex shapes (nanoribbons, rings, user-defined polygons).
- **Modifiers**:
  - Strain fields.
  - Electric and Magnetic fields (Peierls substitution).
  - Disorder (vacancy, onsite energy, hopping).
- **Observables**:
  - Band structures.
  - Local Density of States (LDOS).
  - Transport (Kubo-Greenwood conductivity).
  - Berry curvature (experimentally supported).

## Key Strengths
- **Performance**: The C++ backend ensures that constructing the Hamiltonian and performing KPM calculations is extremely fast, comparable to pure C/Fortran codes.
- **Ease of Use**: The Python API is modern, intuitive, and designed for rapid prototyping (e.g., using decorators for modifiers).
- **Visualization**: Built-in plotting helpers for lattices, bands, and LDOS maps using Matplotlib.

## Inputs & Outputs
- **Inputs**: Python scripts defining lattice vectors, sublattices, and hoppings.
- **Outputs**:
  - `Results` objects containing NumPy arrays of energies, DOS, etc.
  - Plots.

## Interfaces & Ecosystem
- **Dependencies**: NumPy, SciPy, Matplotlib.
- **Interoperability**: Can export matrices to other solvers if needed.

## Performance Characteristics
- **Scaling**: $O(N)$ for KPM methods, allowing systems with $>10^7$ sites on a desktop.
- **Parallelism**: OpenMP multi-threading in the C++ core.

## Comparison with Other Codes
- **vs. Kwant**: Kwant is the standard for *transport* (scattering matrix/Landauer); Pybinding is generally faster for *spectral properties* (DOS/LDOS) of massive systems due to its specialized KPM implementation.
- **vs. TB2J**: TB2J is for magnetism parameters; Pybinding is for electronic structure.

## Application Areas
- **Graphene Nanodevices**: simulating strain and edge effects.
- **Disordered Systems**: Anderson localization studies requiring large statistical ensembles.
- **Moiré Lattices**: Large unit cells of twisted double-bilayer graphene.

## Community and Support
- **Development**: Dean Moldovan.
- **Source**: GitHub.

## Verification & Sources
- **Website**: [https://docs.pybinding.site/](https://docs.pybinding.site/)
- **Verification status**: ✅ VERIFIED
  - Mature and widely used package.
