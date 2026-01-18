# ExaTN

## Official Resources
- **Homepage**: https://github.com/ornl-qci/exatn
- **Repository**: https://github.com/ornl-qci/exatn
- **License**: BSD 3-Clause (Check repository)
- **Institution**: Oak Ridge National Laboratory (ORNL)

## Overview
ExaTN (Exascale Tensor Networks) is a high-performance C++ software library designed for the processing of arbitrary tensor networks on exascale computing platforms. Developed by the Quantum Computing Institute at ORNL, it targets applications that require numerical tensor algebra at the absolute limits of scale, such as quantum many-body theory, coupled cluster methods, and quantum circuit simulation. ExaTN abstracts the complexity of distributed memory and GPU management, providing a unified API for contracting vast tensor networks.

**Scientific domain**: Quantum Chemistry, Many-Body Physics, High-Performance Computing.
**Target user community**: HPC developers, Quantum Chemists, Physics researchers using leadership-class supercomputers.

## Theoretical Methods
- **Tensor Network Contraction**: General algorithm for arbitrary graphs.
- **Tensor Decomposition**: SVD, ISO, etc. on distributed tensors.
- **Optimization**: Graph-based optimization of contraction sequences.
- **Approximation**: Rank-adaptive tensor approximation options.

## Capabilities (CRITICAL)
- **Distributed Computing**: Runs on laptop to leadership clusters (Summit, Frontier).
- **GPU Acceleration**: Native support for NVIDIA (CUDA) and AMD (ROCm) GPUs.
- **APIs**:
    - **Declarative**: High-level graph description.
    - **Executive**: Low-level execution control.
- **Memory Management**: Automatic slicing of large tensors to fit distributed memory.
- **Integration**: Works with `cuQuantum`, `TALSH` (Tensor Algebra Library for Shared Memory).

## Key Strengths
### Exascale Ready
- Designed explicitly for systems with thousands of GPUs.
- Handles tensors larger than the memory of a single node.

### Backend Agnosticism
- Supports multiple runtimes and math libraries (cuTENSOR, rocTENSOR).

## Inputs & Outputs
- **Input**: C++ API calls defining tensor shapes, data, and network connectivity.
- **Output**: Numerical tensor results, contracted scalars.

## Interfaces & Ecosystem
- **Language**: C++ (Primary), Fortran (Bindings exists/planned).
- **Dependencies**: MPI, BLAS, LAPACK, CUDA/ROCm (optional).
- **Integration**: Used as the backend for QCPack (Quantum Circuit Simulation) and other ORNL scientific codes.

## Advanced Features
- **Automated Slicing**: Splits tensors automatically if they don't fit in memory.
- **Async Execution**: Non-blocking tensor operations for latency hiding.

## Performance Characteristics
- **Throughput**: Near-peak FLOP performance on GPUs via vendor libraries.
- **Scaling**: Demonstrated scaling to thousands of nodes on Summit.

## Computational Cost
- **High**: Intended for expensive, large-scale calculations.
- **Efficient**: Maximizes hardware utilization.

## Comparison with Other Codes
- **vs ITensor**: ITensor is more user-friendly for physics/MPS; ExaTN is a lower-level HPC engine for arbitrary massive networks.
- **vs Cyclops (CTF)**: Both are distributed; ExaTN focuses heavily on arbitrary tensor network graphs and GPU backends (cuQuantum).
- **vs TensorNetwork**: TensorNetwork (Google) uses ML backends (TF/JAX); ExaTN uses MPI/C++ backends for traditional HPC environments.

## Application Areas
- **Quantum Circuit Simulation**: Simulating 50+ qubits.
- **Electronic Structure**: Coupled Cluster (CC) theory on tensor networks.
- **Nuclear Physics**: Many-body nuclear structure.

## Best Practices
- **Network Construction**: Define the full network graph before triggering execution to allow global optimization.
- **Resource Allocation**: Careful binding of MPI ranks to GPUs.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/ornl-qci/exatn
2. Paper: "ExaTN: A Scalable Tensor Network Library..." (Frontiers in Applied Math. and Stat.).

**Confidence**: VERIFIED - ORNL project.
**Verification status**: ✅ VERIFIED
