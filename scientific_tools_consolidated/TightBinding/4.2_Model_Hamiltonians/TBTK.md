# TBTK (Tight-Binding ToolKit)

## Official Resources
- **Homepage**: http://www.second-quantization.com/
- **Repository**: https://github.com/dafer45/TBTK
- **License**: Apache-2.0

## Overview
**TBTK** is a C++ library designed for modeling and solving **second-quantized Hamiltonians** on arbitrary discretizable structures. While rooted in tight-binding models, its abstract graph-based architecture allows it to handle a wide variety of quantum mechanical problems, from simple lattices to complex device geometries. It provides a suite of high-performance solvers, including exact diagonalization and the **Kernel Polynomial Method (KPM)**, along with tools for calculating Green's functions and other observables.

**Scientific domain**: Quantum Transport, Second Quantization
**Target user community**: C++ developers simulating mesoscopic devices

## Theoretical Methods
- **Second Quantization**: Abstract representation of creation/annihilation operators on discrete indices.
- **Solvers**:
  - **Diagonalization**: Full spectrum (LAPACK) or partial spectrum (Arnoldi).
  - **Chebyshev (KPM)**: Order-N expansion for Density of States (DOS) and spectral functions.
  - **Block Diagonalization**: Exploiting symmetries.
- **Formalism**: Handles Fermionic and Bosonic statistics (Grand Canonical Ensemble).

## Capabilities
- **Model Construction**:
  - Arbitrary graphs (1D, 2D, 3D, and beyond).
  - Complex indices (subcoordinates, spins, orbitals).
- **Observables**:
  - Density of States (DOS), Local DOS (LDOS).
  - Charge Density, Magnetization.
  - Current Density (Bond currents).
- **Advanced Features**:
  - Self-consistent mean-field loops (Hartree-Fock, Superconductivity).
  - CUDA acceleration for select solvers.

## Key Strengths
- **Generality**: The index system is extremely flexible, allowing models that don't fit into standard "unit cell" descriptions (e.g., quasicrystals, fractals, amorphous systems).
- **Performance**: Written in modern C++ with OpenMP and GPU support, capable of scaling to millions of sites using the Chebyshev solver.
- **Visualization**: Built-in support for exporting data to VTK formats for 3D visualization in Paraview.

## Inputs & Outputs
- **Inputs**: C++ code defining the `Model` object (HoppingAmplitudes).
- **Outputs**:
  - Property containers (Arrays).
  - .vtp/.vtu files for visualization.

## Interfaces & Ecosystem
- **Dependencies**: LAPACK, BLAS, FFTW.
- **Integration**: Designed to be compiled as a shared library.

## Performance Characteristics
- **Efficiency**: State-of-the-art for large sparse systems (KPM).
- **Scalability**: MPI parallelism allows cluster deployment.

## Comparison with Other Codes
- **vs. Kwant**: Kwant is Python-based and focused on scattering. TBTK is C++ and emphasizes the general construction of second-quantized Hamiltonians and spectral properties.
- **vs. Pybinding**: Both use KPM for large systems; TBTK offers a lower-level C++ API which may be preferred for embedding in other high-performance applications.

## Application Areas
- **Mesoscopic Superconductivity**: Modeling Josephson junctions and proximity effects.
- **Quantum Hall Effect**: Edge states in specific geometries.

## Community and Support
- **Development**: Kristofer Björnson.
- **Source**: GitHub.

## Verification & Sources
- **Website**: [http://www.second-quantization.com/](http://www.second-quantization.com/)
- **Verification status**: ✅ VERIFIED
  - Active C++ project.
