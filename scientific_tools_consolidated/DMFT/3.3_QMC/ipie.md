# ipie

## Official Resources
- Homepage: https://github.com/pauxy-qmc/ipie
- Documentation: https://ipie.readthedocs.io/
- Source Repository: https://github.com/pauxy-qmc/ipie
- License: Apache License 2.0

## Overview
**ipie** is a state-of-the-art Python-based Auxiliary-Field Quantum Monte Carlo (AFQMC) package designed for estimating ground state and excited state properties of quantum many-body systems. As a modern successor to PAUXY, it emphasizes performance, modularity, and ease of use. `ipie` is built from the ground up to support high-performance computing on both CPUs and GPUs, making it a powerful tool for *ab initio* quantum chemistry and condensed matter physics.

**Scientific domain**: Quantum Chemistry, Condensed Matter Physics, Electronic Structure
**Target user community**: Researchers in strongly correlated electron systems, quantum chemists needing high-accuracy many-body methods

## Theoretical Methods
- **Auxiliary-Field Quantum Monte Carlo (AFQMC)**: Phaseless AFQMC to control the fermion sign problem.
- **Trial Wavefunctions**: Supports Hartree-Fock and Multi-Slater Determinant (MSD) expansions.
- **Basis Sets**: Compatible with Gaussian-type orbitals (via PySCF) and plane-wave bases.
- **Isomerization Energies**: Capable of chemically accurate energy differences (e.g., < 1 kcal/mol).

## Capabilities (CRITICAL)
- **GPU Acceleration**: Fully optimized NVIDIA GPU support using CuPy and custom CUDA kernels for critical operations.
- **Multi-Slater Determinants**: Efficient GPU-accelerated algorithms for handling MSD trial wavefunctions (10^6+ determinants).
- **Complex Hamiltonians**: Support for complex-valued Hamiltonians relevant for solid-state physics.
- **Automatic Differentiation**: Integration with auto-diff frameworks for property calculations.
- **Scalability**: MPI parallelism for distributed memory and massive parallelization across nodes.

## Key Features

### Performance Optimization:
- **JIT Compilation**: Uses Numba for just-in-time compilation of CPU kernels.
- **CUDA Kernels**: Custom CUDA kernels for MSD operations to maximize throughput on GPUs.

### Quantum Chemistry Interface:
- **PySCF Integration**: Seamlessly generates Hamiltonians and initial guess wavefunctions from PySCF calculations.
- **TrexIO/Dice**: Interfaces with external FCI/SCI codes like Dice and data formats like TrexIO.

### Modular Design:
- **Object-Oriented**: Easy customization of Hamiltonians, propagators, spinors, and estimators.

## Inputs & Outputs
- **Input formats**:
  - Python scripts using the `ipie` API.
  - HDF5 files for Hamiltonian and wavefunction data.
- **Output data types**:
  - Ground state energies and error estimates.
  - One- and two-body reduced density matrices.
  - HDF5 output files for analysis.

## Interfaces & Ecosystem
- **Upstream**: PySCF (integrals/SCF), TrexIO (I/O).
- **Downstream**: Analysis scripts in Python (NumPy/Pandas).
- **acceleration**: CuPy, Numba.

## Workflow and Usage
A typical workflow involves running a mean-field calculation (e.g., in PySCF) to generate integrals and a trial wavefunction. These are processed into `ipie` format. The AFQMC simulation is then launched via a Python driver script, distributing walkers across available MPI ranks or GPUs.

## Performance Characteristics
- **Speed**: Competitive with or faster than C++ codes like QMCPACK/Dice for suitable systems.
- **Accuracy**: Benchmarked to reproduce exact isomerization energies for small basis sets and provide high accuracy for large transition metal complexes.
- **Parallelism**: Hybrid MPI+GPU parallelization.

## Comparison with Other Codes
| Feature | ipie | QMCPACK (AFQMC) | ad_afqmc |
| :--- | :--- | :--- | :--- |
| **Language** | Python (CuPy/Numba) | C++ / CUDA | Python (JAX) |
| **Focus** | Modern, Modular, GPU-first | Production Scale, Integrated | Differentiable Programming |
| **Basis Sets** | Gaussian / Plane-wave | Gaussian / Plane-wave | Gaussian |
| **Performance** | High (GPU Optimized) | High (HPC Optimized) | High (JIT Compiled) |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/pauxy-qmc/ipie
2. "ipie: A Python-Based Auxiliary-Field Quantum Monte Carlo Program" (J. Chem. Theory Comput. article).
3. Documentation: https://ipie.readthedocs.io/

**Verification status**: ✅ VERIFIED
- Source code: OPEN (Apache 2.0)
- Active development: Frequent updates, active maintainers.
- Focus: High-performance Python/GPU AFQMC.
