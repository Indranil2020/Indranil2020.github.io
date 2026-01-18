# TBPLaS

## Official Resources
- **Homepage**: https://www.tbplas.net/
- **Repository**: https://github.com/deepmodeling/tbplas
- **License**: BSD 3-Clause License

## Overview
**TBPLaS** (Tight-Binding Package for Large-scale Simulation) is a high-performance Python package designed for the simulation of **electronic structure** and **quantum transport** in macroscopic tight-binding models. Developed by the **DeepModeling** community, it leverages efficient numerical algorithms—such as the **Tight-Binding Propagation Method (TBPM)** and **Kernel Polynomial Method (KPM)**—to perform calculations on systems with **millions of atomic orbitals**, scaling linearly with system size.

**Scientific domain**: Large-scale Tight-Binding, Quantum Transport, 2D Materials
**Target user community**: Researchers studying disordered systems, Moiré superlattices, and quasicrystals

## Theoretical Methods
- **Tight-Binding Propagation Method (TBPM)**: An $O(N)$ method that calculates time-correlation functions of random states to extract spectral and transport properties without diagonalization.
- **Kernel Polynomial Method (KPM)**: Chebyshev expansion of the density of states and spectral functions.
- **Green's Functions**: Recursive Green's Function (RGF) for exact transport in smaller systems.
- **Exact Diagonalization**: Classic solvers for small systems or band structures.

## Capabilities
- **Observables**:
  - Density of States (DOS) and Local DOS.
  - Optical Conductivity $\sigma(\omega)$.
  - DC Conductivity (Kubo formula).
  - Hall Conductivity ($\sigma_{xy}$) and Chern numbers.
  - Polarization and Dielectric function.
  - Quasieigenstates.
- **Systems**:
  - Graphene, TMDs, and Twisted Bilayers (Moiré).
  - Disordered alloys (Anderson localization).
  - Fractal lattices and quasicrystals.

## Key Strengths
- **Scalability**: Capable of simulating $>10^7$ atoms on a single workstation, thanks to the linear scaling ($O(N)$) of TBPM/KPM.
- **Performance**: Critical kernels are optimized in Cython/Fortran and parallelized with OpenMP/MPI.
- **Ease of Use**: Object-oriented Python API allows for intuitive system construction and analysis.
- **Integration**: Part of the DeepModeling ecosystem (DeepMD-kit), facilitating ML-driven potentials/Hamiltonians.

## Inputs & Outputs
- **Inputs**:
  - Lattice and Hamiltonian definitions (Python objects).
  - Simulation parameters (time steps, energy resolution).
- **Outputs**:
  - Spectral data (DOS, conductivity) as NumPy arrays.
  - Visualization files.

## Interfaces & Ecosystem
- **ASE**: Interface to Atomic Simulation Environment for structure manipulation.
- **DeepModeling**: Potential for future integration with Deep Wannier/Deep Hamiltonian.

## Performance Characteristics
- **Speed**: TBPM is orders of magnitude faster than Exact Diagonalization for large $N$.
- **Efficiency**: Sparse matrix operations reduce memory footprint.

## Comparison with Other Codes
- **vs. KITE**: Both use linear scaling methods (KPM/TBPM). KITE emphasizes "disorder on the fly" and C++ backend; TBPLaS offers a pure Python-centric workflow and includes the powerful propagation method (TBPM).
- **vs. Kwant**: Kwant is better for open systems with leads (scattering); TBPLaS excels at bulk properties of massive disordered systems using spectral methods.

## Application Areas
- **Twistronics**: Electronic properties of twisted bilayer graphene at magic angles (thousands of atoms per cell).
- **Anderson Localization**: Scaling theory of localization in 2D/3D disordered lattices.
- **Quantum Hall Effect**: Calculating Chern numbers in large topological systems.

## Community and Support
- **Development**: DeepModeling Community (Yuan Ping Feng group / Contributors).
- **Source**: GitHub.

## Verification & Sources
- **Website**: [https://www.tbplas.net/](https://www.tbplas.net/)
- **Primary Publication**: Y. Li et al., arXiv:2209.00806 (2022).
- **Verification status**: ✅ VERIFIED
  - Active and modern research tool.
