# SIRIUS

## Official Resources
- Homepage: https://github.com/electronic-structure/SIRIUS
- Documentation: https://sirius-lib.readthedocs.io/
- Source Repository: https://github.com/electronic-structure/SIRIUS
- License: BSD 2-Clause

## Overview
SIRIUS is a domain-specific library for electronic structure calculations, developed at the Swiss National Supercomputing Centre (CSCS). It implements providing a high-performance backend for plane-wave (PP-PW) and full-potential (FP-LAPW) DFT calculations. It is designed from the ground up for hybrid computing architectures (CPU+GPU) and serves as a backend for flagship codes like Quantum ESPRESSO as well as a standalone solver.

**Scientific domain**: High-performance computing, electronic structure
**Target user community**: Developers of DFT codes, HPC centers, researchers needing GPU-accelerated solvers

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets (with Norm-Conserving, Ultrasoft, PAW potentials)
- Full-Potential Linearized Augmented Plane Wave (FP-LAPW)
- Gamma-point specific algorithms
- Iterative solvers (Davidson, Chebyshev)

## Capabilities
- Ground-state electronic structure
- Stress and force calculations
- Magnetization (collinear and non-collinear)
- Spin-orbit coupling
- GPU acceleration (NVIDIA CUDA, AMD ROCm)
- Interface to multiple codes (Quantum ESPRESSO, Elk, Exciting)

## Capabilities
- **Performance**: State-of-the-art GPU utilization.
- **Versatility**: Handles both pseudopotential and all-electron methods in one framework.
- **Interoperability**: Can act as a plugin to accelerate existing Fortran codes.

## Inputs & Outputs
- **Input formats**:
  - JSON-based input files (SIRIUS native)
  - Interface calls from host codes (QE, Elk)
  
- **Output data types**:
  - HDF5 output
  - Standard density/potential files

## Interfaces & Ecosystem
- **Quantum ESPRESSO (QE)**: SIRIUS can replace the internal PW engine of QE.
- **Elk/Exciting**: Provides GPU acceleration for these LAPW codes.
- **CP2K**: Interface in development/testing.

## Performance Characteristics
- **Speed**: Extreme, especially on GPU-dense nodes (e.g., Summit, Piz Daint).
- **Scaling**: Scales to thousands of GPUs. Demonstrates ~160x speedup on NVIDIA Grace Hopper (GH200) vs standard CPU nodes for specific solvers.
- **Memory**: Heavy usage of GPU HBM (High Bandwidth Memory); efficiency depends on keeping data resident on device.

## Computational Cost
- **High Efficiency**: Offloads heavy linear algebra and FFTs to GPU, freeing up CPU for other tasks.
- **Overhead**: Initialization cost is non-trivial; best for large, production-grade systems, not tiny tests.

## Best Practices

### Hybrid Parallelization:
- **Ranks vs Threads**: Use fewer MPI ranks and more OpenMP threads per rank to saturate the GPU.
- **Pools**: Utilize k-point pool parallelization to distribute work across GPUs.

### Installation:
- **Dependencies**: Ensure a functioning CUDA/ROCm toolchain and compatible libraries (SpFFT, SPLA) are installed first.

## Community and Support
- **Hosting**: [GitHub](https://github.com/electronic-structure/SIRIUS).
- **Organization**: Maintained by CSCS (Swiss National Supercomputing Centre).
- **Support**: Active GitHub Issues tracker.

## Verification & Sources
**Primary sources**:
1. Official GitHub: https://github.com/electronic-structure/SIRIUS
2. CSCS Documentation
3. "SIRIUS: domain specific library for electronic structure calculations" (Research citations)

**Confidence**: VERIFIED - Code is a core component of modern HPC electronic structure.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: HPC/Semiconductors/Metals
- Key Feature: GPU Acceleration

## Limitations & Known Constraints
- **Complexity**: As a library, it separates the "physics" from the "solver", requiring the host code to handle high-level workflows.
- **Installation**: Requires a modern C++ toolchain and GPU SDKs.

## Comparison with Other Codes
- **vs Native QE**: SIRIUS is often faster on GPUs due to C++ optimized kernels compared to the native Fortran GPU port of QE.
- **vs BigDFT**: SIRIUS handles standard plane-waves and LAPW; BigDFT handles wavelets.

