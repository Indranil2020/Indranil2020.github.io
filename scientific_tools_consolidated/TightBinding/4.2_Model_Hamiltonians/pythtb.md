# PythTB (Python Tight-Binding)

## Official Resources
- **Homepage**: https://www.physics.rutgers.edu/pythtb/
- **Documentation**: https://www.physics.rutgers.edu/pythtb/usage.html
- **Repository**: https://github.com/pythtb/pythtb
- **License**: MIT License (Assumed, check distribution)

## Overview
**PythTB** is a lightweight Python package for constructing and solving tight-binding models. Developed by David Vanderbilt's group at Rutgers University, it is widely used as a pedagogical tool for teaching **topological band theory** and Berry phase physics. Its simplicity makes it ideal for rapid prototyping of model Hamiltonians, while its specialized routines for calculating **Berry phases**, **Wilson loops**, and **Chern numbers** make it a powerful research tool for topological insulators and semimetals.

**Scientific domain**: Band Theory, Topology, Berry Phase
**Target user community**: Students, Educators, and Researchers in topological matter

## Theoretical Methods
- **Tight-Binding**: Slater-Koster or manual hopping parameters.
- **Topological Invariants**:
  - **Berry Phase**: $\phi = \oint \mathbf{A} \cdot d\mathbf{k}$.
  - **Berry Curvature**: $\Omega(\mathbf{k}) = \nabla \times \mathbf{A}(\mathbf{k})$.
  - **Chern Numbers**: Integration of curvature over the BZ.
  - **$\mathbb{Z}_2$ Invariants**: Via Wilson loop spectra or parity analysis (for inversion symmetric systems).
- **Surface States**: Slab construction to inspect edge modes.

## Capabilities
- **Model Construction**:
  - 0D (molecules) to 3D crystals.
  - Spinful and Spinless fermions.
  - Complex hopping amplitudes.
- **Solvers**:
  - Exact Diagonalization on k-paths or meshes.
- **Analysis**:
  - Band structures.
  - Berry flux through k-space patches.
  - Hybrid Wannier charge centers (Wilson loops).

## Key Strengths
- **Pedagogy**: The syntax is extremely cleaner and intuitive (`model.set_hop(...)`), making it the "Arduino of tight-binding codes."
- **Topology Native**: Unlike general TB codes, PythTB has built-in, robust functions specifically for Berry phases and Wilson loops, reflecting the expertise of the Vanderbilt group.
- **Pure Python**: Zero compilation required; easy to modify and inspect.

## Inputs & Outputs
- **Inputs**: Python scripts.
- **Outputs**:
  - Arrays of eigenvalues/vectors.
  - Matplotlib plots (bands, Berry curvature).

## Interfaces & Ecosystem
- **Dependencies**: NumPy, Matplotlib.
- **Adoption**: Used in widely circulated lecture notes on topological insulators.

## Performance Characteristics
- **Speed**: Pure Python loops can be slow for very large systems or huge k-meshes.
- **Scalability**: Best for small unit cells (model physics). Not designed for large-scale supercells (use Pybinding or Kwant for that).

## Comparison with Other Codes
- **vs. Kwant**: Kwant is for transport/scattering. PythTB is for bulk band topology and Berry phases.
- **vs. Pybinding**: Pybinding is high-performance (C++) for large systems. PythTB is simpler and focused on topological invariants of perfect lattices.

## Application Areas
- **Haldane Model**: The classic example of a Chern insulator.
- **Weyl Semimetals**: Calculating the chirality of Weyl nodes via Berry flux.
- **Education**: Teaching the concept of the geometric phase.

## Community and Support
- **Development**: Rutgers University (Sinisa Coh, David Vanderbilt).
- **Source**: Website / GitHub.

## Verification & Sources
- **Website**: [https://www.physics.rutgers.edu/pythtb/](https://www.physics.rutgers.edu/pythtb/)
- **Verification status**: ✅ VERIFIED
  - The standard pedagogical tool for topological bands.
