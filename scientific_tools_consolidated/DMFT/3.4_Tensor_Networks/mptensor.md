# mptensor

## Official Resources
- **Repository**: https://github.com/smorita/mptensor
- **License**: GNU General Public License v3.0
- **Developer**: Satoshi Morita (ISSP, University of Tokyo)

## Overview
mptensor is a parallel C++ library for tensor calculations, designed to provide a high-performance backend for tensor network simulations on distributed memory systems. It mimics the user-friendly interface of NumPy/SciPy while executing operations in parallel across thousands of cores using MPI and OpenMP. mptensor serves as the core engine for the `TeNeS` solver, enabling it to scale to leadership-class supercomputers.

**Scientific domain**: Tensor Networks, High-Performance Computing.
**Target user community**: Developers of C++ tensor network codes (e.g., DMRG, PEPS) requiring MPI parallelism.

## Theoretical Methods
- **Tensor Algebra**: Contraction, SVD, QR, Reshape, Transpose.
- **Distributed Computing**: Block-cyclic distribution of tensors over MPI processes.
- **Symmetries**: U(1) and Z2 symmetry support (block-sparse arithmetic).

## Capabilities (CRITICAL)
- **NumPy-like C++ API**: Intuitive syntax for C++ developers familiar with Python.
- **Massive Parallelism**: Hybrid MPI/OpenMP parallelization.
- **ScaLAPACK Integration**: Uses ScaLAPACK for distributed linear algebra (SVD, Eigen).
- **Routines**: Full suite of tensor contraction (`tensordot`), decompositions, and element-wise operations.
- **Micro-optimization**: efficient local tensor operations.

## Key Strengths
### Ease of Use + Speed
- Writes like Python, runs like optimized HPC code.
- Abstracts away the complexity of MPI communication for tensor slicing.

### Scalability
- Proven scaling on massive partitions of the Fugaku supercomputer (via TeNeS).

## Inputs & Outputs
- **Input**: C++ source code linking the library.
- **Output**: Distributed tensor objects, binaries.

## Interfaces & Ecosystem
- **Dependencies**: MPI, LAPACK, ScaLAPACK.
- **Integration**: Used by `TeNeS`, `Tensordot`.
- **Language**: C++11 (Template library).

## Advanced Features
- **Symmetry Sectors**: Automatic handling of quantum number conservation.
- **Tensordot Code Geeration**: Helper tool to generate optimized contraction kernels.

## Performance Characteristics
- **Efficiency**: Near-optimal usage of distributed linear algebra libraries.
- **Overhead**: Minimal communication overhead for large block sizes.

## Computational Cost
- **High**: Designed for heavy-lifting tasks; lightweight compared to the physics simulations it enables.

## Comparison with Other Codes
- **vs CTF (Cyclops)**: Both are distributed; mptensor focuses on a NumPy-like API for specific physics use cases (TeNeS).
- **vs ITensor**: ITensor (v3) uses OpenMP but distributed MPI is less central/native than in mptensor's design philosophy.

## Application Areas
- **Tensor Network Solvers**: Backend for 2D/3D TN algorithms.
- **RG Methods**: Tensor Renormalization Group (TRG, HOTRG).

## Best Practices
- **Linking**: Ensure optimized vendor BLAS/ScaLAPACK (Intel MKL, Cray LibSci) are used.
- **Process Grid**: Match MPI process grid to physical hardware topology for best communication.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/smorita/mptensor
2. University of Tokyo ISSP activity reports.

**Confidence**: VERIFIED - Active and used in production codes.
**Verification status**: ✅ VERIFIED
