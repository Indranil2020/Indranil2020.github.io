# EDRIXS

## Official Resources
- **Homepage**: https://github.com/EDRIXS/edrixs
- **Documentation**: https://nsls-ii.github.io/edrixs/
- **Source Repository**: https://github.com/EDRIXS/edrixs
- **PyPI**: https://pypi.org/project/edrixs/
- **License**: BSD-3-Clause

## Overview
**EDRIXS** (Exact Diagonalization for Resonant Inelastic X-ray Scattering) is an open-source toolkit designed for simulating X-ray Absorption Spectroscopy (XAS), Resonant Inelastic X-ray Scattering (RIXS), and Resonant Magnetic X-ray Scattering (RMXS). Built on exact diagonalization (ED) of model Hamiltonians, it is particularly suited for studying strongly correlated materials where local interactions play a critical role. EDRIXS combines a high-performance Fortran 90 core for heavy numerical tasks with a flexible Python interface for setup and post-processing, allowing users to handle complex interactions like spin-orbit coupling, crystal structure fields, and Coulomb interactions within a user-defined cluster.

**Scientific domain**: Strongly correlated electron systems, X-ray spectroscopy (XAS, RIXS, RMXS), Transition metal oxides, f-electron systems
**Target user community**: Researchers in condensed matter physics and physical chemistry focusing on spectroscopic analysis of correlated materials

## Theoretical Methods
- **Exact Diagonalization (ED)**: Solves many-body Hamiltonians for eigenvalues and eigenvectors
- **Lanczos Algorithm**: Efficient iterative method for finding ground and excited states
- **Krylov Subspace Techniques**: Used for calculating spectral functions and response functions
- **Cluster Perturbation Theory**: For treating extended systems beyond small clusters
- **Model Hamiltonians**:
  - Anderson Impurity Model
  - Hubbard Model parameters
  - Crystal Field Theory
  - Spin-Orbit Coupling
  - Slater-Condon parameters for Coulomb interactions
- **Wannier Functions**: Integration with Wannier90 to define realistic model parameters from DFT

## Capabilities
- **Spectroscopy Simulation**:
  - X-ray Absorption Spectroscopy (XAS)
  - Resonant Inelastic X-ray Scattering (RIXS)
  - Resonant Magnetic X-ray Scattering (RMXS)
- **Hamiltonian Solvers**:
  - **Fortran Solver**: MPI-parallelized for large Hilbert spaces (based on ARPACK)
  - **Python Solver**: Pure Python implementation for small systems (Hilbert space < ~10,000)
- **Complex Interaction Modeling**:
  - Full multiplet theory
  - Charge transfer effects
  - Low-symmetry crystal fields
- **Transition Amplitudes**: Calculation of non-spin-flip and spin-flip processes

## Key Strengths
- **Hybrid Architecture**: Combines the ease of use of Python with the performance of Fortran/MPI.
- ** versatility**: Applicable to single atoms, small clusters, and impurity models.
- **Spectroscopy Focus**: Specialized for simulating modern resonant X-ray experiments.
- **Integration**: Seamlessly utilizes parameters from DFT+Wannier90 or DMFT calculations.
- **Open Source**: Community-driven development with BSD licensing.

## Inputs & Outputs
- **Inputs**:
  - Electronic structure parameters (hopping integrals, Coulomb U, J)
  - Geometry and cluster definitions
  - Incident/outgoing photon energies and polarizations
  - Slater-Condon parameters
- **Outputs**:
  - Spectral functions (Intensity vs Energy)
  - Eigenvalues and Eigenvectors
  - Expectation values of operators
  - Transition amplitudes

## Interfaces & Ecosystem
- **Python Interface**: comprehensive API for workflow management, pre-processing, and plotting (`import edrixs`).
- **Wannier90**: Can read and use Wannier functions to construct realistic tight-binding Hamiltonians.
- **DFT Integration**: Bridges first-principles calculations (e.g., via Wannier90) with many-body model simulations.
- **DMFT**: Compatible with Dynamical Mean-Field Theory workflows for impurity problems.

## Performance Characteristics
- **Parallelization**: MPI-based parallelism for the Fortran solver allows scaling to multi-core clusters.
- **Efficiency**: Optimized Lanczos and Krylov solvers for sparse matrices.
- **Scalability**: Limited by the exponential scaling of the Hilbert space size typical of Exact Diagonalization methods.
- **Memory**: Efficient handling of sparse Hamiltonian matrices.

## Limitations & Known Constraints
- **System Size**: Restricted to small clusters or impurity models due to exponential scaling of ED.
- **Python Solver**: "Pure" Python solver is limited to small Hilbert spaces (< 10,000 basis states); larger problems require the Fortran/MPI backend.
- **Approximation**: Cluster methods may miss long-range correlations present in bulk systems (though CPT helps).

## Comparison with Other Codes
- **vs. Quanty**: Both perform multiplet calculations; EDRIXS is open-source (BSD) and emphasizes Python/Wannier integration.
- **vs. ALPS**: ALPS is a general framework for lattice models; EDRIXS is specifically optimized for X-ray spectroscopy workflows (XAS/RIXS).
- **vs. CTM4XAS**: CTM4XAS is often GUI-based; EDRIXS offers a programmable Python environment for advanced users.

## Application Areas
- **Transition Metal Oxides**: Studying d-orbital physics, magnetism, and charge transfer.
- **Lanthanides/Actinides**: Modeling f-electron systems with strong spin-orbit coupling.
- **High-Tc Superconductors**: Analyzing RIXS spectra to understand magnetic and orbital excitations.
- **Quantum Materials**: Investigating topological and correlated phases via spectroscopic signatures.

## Community and Support
- **Source**: Hosted on GitHub with issue tracking and contributions.
- **Documentation**: Hosted on GitHub Pages (NSLS-II).
- **Development**: Originally developed at Brookhaven National Laboratory (COMSCOPE project).

## Verification & Sources
- **Official Website**: [https://github.com/EDRIXS/edrixs](https://github.com/EDRIXS/edrixs)
- **Documentation**: [https://nsls-ii.github.io/edrixs/](https://nsls-ii.github.io/edrixs/)
- **Primary Publication**: Wang, Y. L., et. al., "EDRIXS: An open source toolkit for simulating RIXS spectra...". Comput. Phys. Commun. 243, 151 (2019).
- **Verification status**: ✅ VERIFIED
  - Code is active and open source.
  - Validated against standard atomic multiplet codes and experimental data.
