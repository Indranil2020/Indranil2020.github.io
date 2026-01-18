# Quantica.jl

## Official Resources
- **Homepage**: https://pablosanjose.github.io/Quantica.jl/stable/
- **Repository**: https://github.com/pablosanjose/Quantica.jl
- **License**: MIT License

## Overview
**Quantica.jl** is a high-performance **Julia** framework for the construction and simulation of quantum lattice systems. Designed as a modern, faster alternative to Python-based tools like Kwant, it provides an expressive API for defining tight-binding Hamiltonians and efficiently calculating spectral and transport properties using **Green's function** methods. It natively supports **superconducting systems** (Bogoliubov-de Gennes) and allows for arbitrary parametric dependence of Hamiltonians.

**Scientific domain**: Mesoscopic Physics, Topological Superconductivity, Quantum Transport
**Target user community**: Theorists aiming for high-performance simulations of tight-binding models

## Theoretical Methods
- **Tight-Binding & BdG**: Supports standard tight-binding and superconducting Hamiltonians in Nambu space.
- **Recursive Green's Function (RGF)**: Efficiently computes transport (S-matrix) and local properties for quasi-1D systems (leads + scattering region).
- **Kernel Polynomial Method (KPM)**: (Via extensions or integration) for spectral properties of large systems.
- **Sparse Diagonalization**: Fast eigenvalue solvers for band structures.

## Capabilities
- **System Building**:
  - "Builder" pattern (similar to Kwant) for defining lattices, hoppings, and shapes.
  - Parametric Hamiltonians ($H(t, B, \dots)$) without recompilation.
- **Observables**:
  - Local Density of States (LDOS).
  - Josephson Currents ($I(\phi)$).
  - Transmission and Conductance.
  - Band structures.
- **Physics**:
  - Majorana fermions in nanowires.
  - Quantum spin Hall effect.
  - Andreev reflection.

## Key Strengths
- **Performance**: Written in pure Julia, it benefits from JIT compilation, often outperforming mixed Python/C codes for Hamiltonian generation and custom loops.
- **Superconductivity**: First-class support for Nambu spinors and BdG physics, simplifying the study of hybrid superconductor-semiconductor devices.
- **Expressiveness**: Concise, mathematical syntax for defining models.

## Inputs & Outputs
- **Inputs**: Julia scripts using the `Quantica` DSL.
- **Outputs**:
  - Julia structs (Green's functions).
  - Plotting recipes for `Makie.jl` or `Plots.jl`.

## Interfaces & Ecosystem
- **Julia Ecosystem**: Interoperable with `LinearAlgebra`, `SparseArrays`, `KrylovKit` (diagonalization).
- **Visualisation**: Native plotting recipes for visualizing lattices and fields.

## Performance Characteristics
- **Speed**: Hamiltonian construction is extremely fast. RGF solver is comparable to optimized Fortran/C codes.
- **Scalability**: Capable of handling systems with $10^5-10^6$ orbitals on a single node.

## Comparison with Other Codes
- **vs. Kwant**: Quantica is the "Julia answer" to Kwant. It is faster for constructing Hamiltonians and iterating over parameters, but Kwant has a mature, larger ecosystem (Tkwant, etc.).
- **vs. PyBinding**: Quantica offers more advanced transport capabilities (Green's functions) beyond just band structure.

## Application Areas
- **Topological Quantum Computing**: Modeling Majorana zero modes in superconductor-semiconductor heterostructures.
- **Josephson Junctions**: Current-phase relationships in complex geometries.
- **Twisted Bilayers**: Moiré Hamiltonians (performance benefit for large unit cells).

## Community and Support
- **Development**: Pablo San-Jose (ICMM-CSIC, Madrid).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/pablosanjose/Quantica.jl](https://github.com/pablosanjose/Quantica.jl)
- **Verification status**: ✅ VERIFIED
  - Active and highly regarded in the Julia physics community.
