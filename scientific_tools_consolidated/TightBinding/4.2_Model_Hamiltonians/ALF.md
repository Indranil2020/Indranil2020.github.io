# ALF

## Official Resources
- **Homepage**: https://alf.physik.uni-wuerzburg.de/
- **Repository**: https://git.physik.uni-wuerzburg.de/ALF/ALF
- **Documentation**: https://alf.physik.uni-wuerzburg.de/Documentation/
- **License**: GPL-3.0

## Overview
**ALF** (Algorithms for Lattice Fermions) is a high-performance, open-source software package designed for **Quantum Monte Carlo (QMC)** simulations of strongly correlated fermion systems. It specifically implements the **Auxiliary-Field QMC (AFQMC)** method (also known as Determinantal QMC) to solve a wide class of lattice models at finite temperature or in the ground state (projective). Its key strength lies in its generality: users can define arbitrary Hamiltonians writable in terms of single-body and squared single-body operators.

**Scientific domain**: Strongly Correlated Electrons, Lattice Gauge Theory
**Target user community**: Theorists studying Hubbard, Kondo, and topological phase transitions

## Theoretical Methods
- **Method**: Auxiliary-Field Quantum Monte Carlo (AFQMC) / Determinantal QMC.
- **Transformations**: Trotter-Suzuki decomposition for imaginary time and Hubbard-Stratonovich transformation to decouple interactions.
- **Ensembles**:
  - Finite Temperature (Grand Canonical).
  - Projective (Zero Temperature / Canonical).
- **Stabilization**: Matrix decomposition stabilization for long imaginary time propagation.

## Capabilities
- **Models**:
  - Hubbard Model (single/multi-orbital).
  - Kondo Lattice Model.
  - Periodic Anderson Model.
  - $Z_2$ Lattice Gauge Theories.
  - User-defined generic Hamiltonians.
- **Observables**:
  - Equal-time and time-displaced Green's functions.
  - Spin and charge susceptibilities.
  - Renyi Entanglement Entropy.
  - Superconducting pairing correlations.
- **Analysis**: Stochastic Maximum Entropy (MaxEnt) for analytic continuation to real frequencies.

## Key Strengths
- **Versatility**: Unlike many QMC codes hardwired for the Hubbard model, ALF provides a generic interface (`Hamiltonian_interface`) to defining any suitable interacting model.
- **Optimized**: Comparison-free Fortran implementation with MPI parallelism for massive sampling.
- **Documentation**: Extensive documentation and tutorials provided by the Würzburg group.

## Inputs & Outputs
- **Inputs**:
  - `parameters` file: Lattice size, temperature, interaction strength $U$.
  - Hamiltonian definition (Fortran modules).
- **Outputs**:
  - Data files for observables (bins for error analysis).
  - Python scripts (`pyALF`) for post-processing and plotting.

## Interfaces & Ecosystem
- **pyALF**: A Python interface to manage simulation workflows, generating input files and analyzing results.
- **HDF5**: Optional support for structured data output.

## Performance Characteristics
- **Scaling**:
  - Linear with inverse temperature $\beta$.
  - Cubic $O(N^3)$ with system volume $N$ (number of orbitals).
- **Parallelism**: MPI parallelization over Markov chains (embarrassingly parallel sampling).

## Comparison with Other Codes
- **vs. QUEST**: QUEST is another major DQMC code. ALF is arguably more modern in its software engineering (Fortran 2003 objects) and offers a more generalized Hamiltonian interface.
- **vs. ALPS**: ALPS offers a suite of solvers (including QMC); ALF is a specialized, deep tool for AFQMC with features like projective algorithms and entanglement entropy that standard ALPS applications might lack.

## Application Areas
- **Quantum Criticality**: Studying phase transitions in heavy fermion systems.
- **Topological Phases**: Simulating topological insulators with interactions (Kane-Mele-Hubbard).
- **Graphene**: Hubbard model on honeycomb lattices.

## Community and Support
- **Development**: University of Würzburg (Fakher Assaad group).
- **Source**: GitLab (University hosted).

## Verification & Sources
- **Website**: [https://alf.physik.uni-wuerzburg.de/](https://alf.physik.uni-wuerzburg.de/)
- **Primary Publication**: T. C. Lang et al., "ALF: The algorithms for lattice fermions package", SciPost Phys. (2019).
- **Verification status**: ✅ VERIFIED
  - Active and well-supported community code.
