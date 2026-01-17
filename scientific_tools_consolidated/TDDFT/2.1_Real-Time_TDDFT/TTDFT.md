# TTDFT

## Official Resources
- **Repository**: https://github.com/ttdftdev/ttdft_public
- **License**: Information available in repository
- **Documentation**: See GitHub README

## Overview
TTDFT is a high-performance Real-Space Time-Dependent Density Functional Theory code designed for modern heterogeneous computing architectures. It is written in C++ and leverages NVIDIA GPUs for significant acceleration of key operations, such as matrix multiplications in the Kohn-Sham propagation.

**Scientific domain**: Molecular dynamics, electronic structure, excited states
**Target user community**: HPC users with GPU resources, developers of real-space methods

## Theoretical Methods
- **Real-Space TDDFT**: Grid-based discretization of the Kohn-Sham equations
- **Tensor Structured Algorithms**: Utilization of Tucker tensor decomposition
- **Tamm-Dancoff Approximation (TDA)**: Supported
- **Real-Time Propagation**: Magnus expansion and other integrators
- **Chebyshev Filtering**: Used for efficient eigensolvers

## Capabilities
- **GPU Acceleration**: Native support for CUDA, cuBLAS, cuSparse
- **Real-Time Dynamics**: Simulation of electron dynamics under external fields
- **Large-Scale Systems**: Optimized for systems that benefit from tensor compression
- **Ground State**: Fast ground state convergence via Chebyshev filtering

## Inputs & Outputs
- **Input formats**: text-based input files (implied from standard C++ scientific codes)
- **Output data types**:
  - Energy and forces
  - Time-dependent dipole moments
  - Absorption spectra (via Fourier transform of dipoles)
  - Density cubes

## Performance Characteristics
- **Speed**: ~8x speedup reported for matrix-matrix multiplications on GPUs compared to CPU-only.
- **Parallelization**: Hybrid MPI and CUDA.
- **Efficiency**: Tensor-structured operations reduce memory footprint and computational complexity.

## Computational Cost
- **Memory**: Tensor decomposition helps reduce storage requirements for large grids.
- **Hardware**: Requires NVIDIA GPUs for optimal performance (module load cuda).

## Limitations & Known Constraints
- **Hardware**: Strongly tied to NVIDIA ecosystem (CUDA).
- **Compilation**: Heterogeneous build process requires careful environment setup (MPI + NVCC).
- **Documentation**: Primary documentation is the GitHub README and associated papers.

## Comparison with Other Codes
- **vs Octopus**: Both real-space; TTDFT emphasizes tensor compression and C++/CUDA acceleration.
- **vs GCEED**: GCEED focuses on coupled Maxwell-TDDFT; TTDFT focuses on tensor algorithms and GPU compute.

## Best Practices
- **Compilation**: Ensure `cuda` modules are loaded.
- **GPU Usage**: Target problems large enough to saturate GPU compute capability for best efficiency.

## Citations
- **Primary**: "TTDFT: A GPU accelerated Tucker tensor DFT code for large-scale Kohn-Sham DFT calculations" (arXiv/OSTI).
