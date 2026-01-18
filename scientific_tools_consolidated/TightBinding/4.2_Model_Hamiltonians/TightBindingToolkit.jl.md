# TightBindingToolkit.jl

## Official Resources
- **Homepage**: https://anjishnubose.github.io/TightBindingToolkit.jl/
- **Repository**: https://github.com/Anjishnubose/TightBindingToolkit.jl
- **License**: MIT License

## Overview
**TightBindingToolkit.jl** is a feature-rich Julia package for the construction, solution, and analysis of generic tight-binding models. It excels in the study of **topological phases of matter**, providing built-in tools for Berry curvature, Chern numbers, and Majorana modes in superconductors. It supports both standard electronic Hamiltonians and Bogoliubov-de Gennes (BdG) Hamiltonians for superconductors, making it a versatile tool for defining custom lattice models in 1D, 2D, and 3D.

**Scientific domain**: Topological Insulators/Superconductors, Band Theory
**Target user community**: Theorists exploring topological phases and lattice models

## Theoretical Methods
- **Tight-Binding**: General hopping on Bravais lattices.
- **BdG Formalism**: Particle-hole symmetric Hamiltonians for superconductivity ($H_{BdG} = \begin{pmatrix} H_0 & \Delta \\ \Delta^\dagger & -H_0^* \end{pmatrix}$).
- **Topology**: 
  - Berry connection $\mathbf{A}(\mathbf{k}) = -i \langle u | \nabla_k | u \rangle$.
  - Berry curvature and Chern numbers.
- **Green's Functions**: Momentum-space Green's functions $G(\omega, \mathbf{k})$.

## Capabilities
- **Model Construction**:
  - Arbitrary unit cells and hopping ranges.
  - Multi-orbital bases.
- **Calculations**:
  - Band structures and DOS.
  - 2D Fermi surfaces / Constant energy contours.
  - Topological invariants (Chern number, winding number).
  - Magnetic susceptibility $\chi(\mathbf{q})$.
- **Advanced**:
  - Flux insertion (Peierls substitution) for Hofstadter butterfly spectra.

## Key Strengths
- **Topological Toolkit**: Unlike generic TB codes, it has specific high-level functions for topological invariants, saving the user from implementing Berry phase integration manually.
- **BdG Support**: Native handling of superconducting pairing terms, essential for studying topological superconductors and Majorana fermions.
- **Julia Efficiency**: Fast numerical diagonalization and integration, suitable for parameter sweeps phase diagrams.

## Inputs & Outputs
- **Inputs**: Julia scripts defining lattice vectors, orbitals, and hoppings.
- **Outputs**:
  - Plot objects (Plots.jl/Makie).
  - Data arrays.

## Interfaces & Ecosystem
- **Dependencies**: `LinearAlgebra`, `Combinatorics`.
- **Plotting**: Integrated recipes for standard Julia plotting libraries.

## Performance Characteristics
- **Speed**: High performance for dense k-grids due to Julia's compilation.
- **Parallelism**: Threaded loops for k-space integration.

## Comparison with Other Codes
- **vs. PythTB**: Similar scope, but TightBindingToolkit.jl leverages Julia's speed and has deeper support for BdG/superconductivity.
- **vs. Quantica.jl**: Quantica is another Julia TB code; TightBindingToolkit.jl is perhaps more focused on the *analysis* (susceptibility, topology) of bulk Hamiltonians rather than device transport.

## Application Areas
- **Topological Superconductivity**: Searching for Majorana zero modes in nanowires.
- **Quantum Anomalous Hall**: Studying Chern insulators on honeycomb lattices.
- **Susceptibility**: Calculating nesting vectors in Fermi surfaces.

## Community and Support
- **Development**: Anjishnu Bose.
- **Source**: GitHub.

## Verification & Sources
- **Documentation**: [https://anjishnubose.github.io/TightBindingToolkit.jl/dev/](https://anjishnubose.github.io/TightBindingToolkit.jl/dev/)
- **Verification status**: ✅ VERIFIED
  - Active research package.
