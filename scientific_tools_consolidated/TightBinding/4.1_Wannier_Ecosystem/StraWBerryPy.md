# StraWBerryPy

## Official Resources
- **Repository**: https://github.com/strawberrypy-developers/strawberrypy
- **Documentation**: https://strawberrypy.readthedocs.io/
- **License**: MIT License

## Overview
**StraWBerryPy** (Single-poinT and local invaRiAnts for Wannier Berriologies in Python) is a specialized Python package for calculating topological invariants and quantum geometrical properties in **non-crystalline** and **disordered** topological insulators. Unlike standard tools that rely on Bloch's theorem and Brillouin zone integration, StraWBerryPy operates in real space, making it uniquely capable of characterizing topological phases in amorphous materials, quasicrystals, and systems with broken translational symmetry.

**Scientific domain**: Topological Insulators, Disordered Systems, Amorphous Materials, Quantum Geometry
**Target user community**: Researchers studying topology in non-crystalline matter and disordered solids

## Theoretical Methods
- **Real-Space Topology**: Approaches that do not require k-space integration.
- **Single-Point Formulas**: Efficient calculation of invariants using data from a single point in the supercell configuration space.
- **Local Markers**:
  - **Local Chern Marker (LCM)**: For Time-Reversal Symmetry broken systems (e.g., Quantum Anomalous Hall).
  - **Local Spin Chern Marker (LSCM)**: For Time-Reversal Invariant systems (e.g., Quantum Spin Hall).
  - **Localization Marker**: To distinguish trivial from topological insulating states locally.
- **Supercell Framework**: Handles periodic boundary conditions with large unit cells to simulate disorder.

## Capabilities
- **Topological Invariant Calculation**:
  - Chern Number ($C$)
  - Spin Chern Number ($C_s$)
  - $\mathbb{Z}_2$ invariants (via spin Chern numbers)
- **Local Analysis**: Spatially resolved calculation of topological markers to identify phase boundaries or defects.
- **Disorder Handling**: Built-in tools to add on-site or hopping disorder to clean models.
- **Boundary Conditions**: Supports both Periodic Boundary Conditions (PBC) for bulk invariants and Open Boundary Conditions (OBC) for finite clusters.

## Key Strengths
- **Non-Crystalline Support**: One of the few tools dedicated to topology without translational symmetry.
- **Computational Efficiency**: Single-point formulas are often cheaper than full Brillouin zone integration for large supercells.
- **Interface Flexibility**: Works with multiple tight-binding backends (PythTB, TBmodels).
- **Ab Initio Readiness**: Can analyze realistic materials via Wannier90 Hamiltonians (through WannierBerri).

## Inputs & Outputs
- **Inputs**:
  - Tight-binding models (PythTB/TBmodels objects)
  - Wannier90 Hamiltonian files (`*_hr.dat`)
  - Disorder parameters (strength, type)
- **Outputs**:
  - Global topological invariants (Integer values)
  - Local marker maps (Spatial distribution of topology)
  - Quantum geometric tensor components

## Interfaces & Ecosystem
- **PythTB**: Native support for models built with the Python Tight-Binding package.
- **TBmodels**: Interface for reading/writing models in TBmodels format.
- **WannierBerri**: Integration for reading *ab initio* Wannier Hamiltonians.
- **Wannier90**: Indirect support via WannierBerri for realistic materials analysis.

## Performance Characteristics
- **Scaling**: Efficient for large supercells due to real-space formulation.
- **Sparse Algebra**: Utilizes sparse matrix operations for Hamiltonian diagonalization.
- **Parameter Sweeps**: optimized for calculating phase diagrams across disorder strengths.

## Limitations & Known Constraints
- **Supercell Size**: Accuracy of single-point formulas improves with supercell size; small cells may have finite-size effects.
- **3D Topology**: Primarily focused on 2D invariants (Chern/Spin-Chern) and their dimensional reductions; 3D strong indices may require specific setups.

## Comparison with Other Codes
- **vs. Z2Pack**: Z2Pack uses Wilson loops in k-space (requires periodicity); StraWBerryPy works in real-space (no periodicity required).
- **vs. WannierTools**: WannierTools computes surface states and invariants mostly in k-space (slab method); StraWBerryPy offers local real-space markers.
- **vs. PythTB**: PythTB provides the model construction; StraWBerryPy adds the advanced topological analysis layer.

## Application Areas
- **Amorphous Topological Insulators**: Materials with topological properties but no long-range order.
- **Disordered Systems**: Studying the robustness of edge states against impurities (Anderson localization).
- **Quasicrystals**: Topological phases in quasi-periodic lattices.
- **Defects**: Analyzing local topology around vacancies or domain walls.

## Community and Support
- **Source**: Open development on GitHub.
- **Documentation**: Tutorials available on ReadTheDocs.
- **Issues**: Tracker for bug reports and feature requests.

## Verification & Sources
- **Repository**: [https://github.com/strawberrypy-developers/strawberrypy](https://github.com/strawberrypy-developers/strawberrypy)
- **Documentation**: [https://strawberrypy.readthedocs.io](https://strawberrypy.readthedocs.io)
- **Verification status**: ✅ VERIFIED
  - Active codebase.
  - Methodology based on established real-space topology markers (Kitaev, Prodan, etc.).
