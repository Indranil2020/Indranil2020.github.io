# TightBinding++ (TB++)

## Official Resources
- **Homepage**: https://tightbinding.github.io/
- **Repository**: https://github.com/TightBinding/tbpp
- **License**: MPL-2.0

## Overview
**TightBinding++** is a modern C++11 framework for efficient simulation of tight-binding models. Designed with modularity and performance in mind, it provides a unified interface to compute a wide range of electronic properties, from standard band structures to advanced topological invariants like Berry curvature and Chern numbers. It supports parallel execution via OpenMP and HDF5 I/O, making it suitable for both rapid prototyping and data-intensive research.

**Scientific domain**: Band Theory, Topological Matter
**Target user community**: C++ developers needing a TB library backend

## Theoretical Methods
- **Tight-Binding**: General orthogonal tight-binding models.
- **Linear Response**: Kubo formula for optical conductivity $\sigma_{\alpha\beta}(\omega)$.
- **Topology**:
  - Berry connection and curvature $\Omega(\mathbf{k})$.
  - Chern numbers via integration over the Brillouin Zone.
- **Green's Functions**: Iterative methods for DOS and transport.

## Capabilities
- **Model Construction**:
  - Simple API for adding sites and hoppings.
  - 1D, 2D, 3D lattices.
- **Observables**:
  - Band structures (complex and real).
  - Density of States (DOS).
  - Berry phase properties.
- **I/O**:
  - HDF5 support for storing large result sets.
  - Python bindings (experimental) allow scripting.

## Key Strengths
- **Performance**: Written in highly optimized C++11 with OpenMP, ensuring fast execution on multi-core machines.
- **Modern Design**: Uses modern C++ paradigms (RAII, smart pointers) making the library safer and easier to use than legacy Fortran/C codes.
- **Topology Native**: Unlike many older TB codes, topological quantities (Berry curvature) are first-class citizens in the API.

## Inputs & Outputs
- **Inputs**: C++ driver code or input parameter files.
- **Outputs**: HDF5 files containing bands, curvature, etc.

## Interfaces & Ecosystem
- **Dependencies**: Eigen3, HDF5, FFTW.
- **Visualization**: Helper scripts for plotting results.

## Performance Characteristics
- **Efficiency**: Matrix operations are vectorized via Eigen.
- **Parallelism**: OpenMP threading for k-point loops.

## Comparison with Other Codes
- **vs. TBTK**: Both are C++ TB libraries. TightBinding++ has a slightly stronger focus on standard band/topology tasks, while TBTK focuses on general graph-based models and Green's functions.
- **vs. PythTB**: PythTB is Python (slower, easier); TightBinding++ is C++ (faster, steeper learning curve).

## Application Areas
- **Topological Photovoltaics**: Calculating shift currents (via Kubo).
- **Material Screening**: Fast band structure generation for databases.

## Community and Support
- **Development**: TightBinding team (GitHub).
- **Source**: GitHub.

## Verification & Sources
- **Website**: [https://tightbinding.github.io/](https://tightbinding.github.io/)
- **Verification status**: ✅ VERIFIED
  - Valid C++ framework.
