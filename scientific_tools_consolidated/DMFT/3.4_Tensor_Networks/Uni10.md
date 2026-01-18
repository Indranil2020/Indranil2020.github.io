# Uni10

## Official Resources
- **Homepage**: https://uni10.gitlab.io/
- **Repository**: https://github.com/yingjerkao/uni10
- **License**: LGPL (Lesser GNU Public License)
- **Developers**: Ying-Jer Kao Group (National Taiwan University)

## Overview
Uni10 (Universal Tensor Network Library) is an open-source C++ library designed to facilitate the development and implementation of tensor network algorithms. It balances high-performance C++ execution with ease of use, providing a "Network" class to manage complex tensor diagrams intuitively. Uni10 is particularly focused on enabling researchers to write readable code for algorithms like DMRG, PEPS, and MERA while benefiting from behind-the-scenes optimization and GPU acceleration.

**Scientific domain**: Quantum Many-Body Physics, Tensor Networks.
**Target user community**: Researchers developing custom tensor network codes (C++ or Python).

## Theoretical Methods
- **Tensor Operations**: Contraction, Permutation, Reshaping, Bond fusion.
- **Decompositions**: SVD, QR, Eigendecomposition, LQ.
- **Algorithmic Support**: Primitives for building DMRG, iTEBD, PEPS, MERA.
- **Symmetries**: Support for Abelian symmetries (U(1), Z2) via block-sparse tensors.

## Capabilities (CRITICAL)
- **Network Class**: Object-oriented management of tensor networks; handles index bookkeeping automatically.
- **Python Wrappers**: `pyUni10` allows rapid prototyping in Python with C++ speed.
- **GPU Support**: (Uni10 v2+) Support for CUDA acceleration.
- **Fermions**: Handling of fermionic swap gates and statistics.
- **Optimized Contraction**: Heuristic algorithms to determine optimal pairwise contraction orders.

## Key Strengths
### Usability vs Performance
- Bridges the gap between high-level interpreted code (Python) and bare-metal C++.
- Encapsulates index management, reducing manual bookkeeping errors.

### 2D Focus
- Specifically designed with PEPS and MERA in mind, capable of handling higher-order tensors.

## Inputs & Outputs
- **Input**: C++/Python scripts constructing tensors and networks.
- **Output**: Contracted values, Optimized tensors, Thermodynamic properties.

## Interfaces & Ecosystem
- **Language**: C++ (Core), Python (Bindings).
- **Dependencies**: BLAS, LAPACK.
- **Documentation**: API references and tutorials for DMRG/MERA.

## Advanced Features
- **Block-Sparse Tensors**: Efficient storage and computation for symmetric systems.
- **Network Optimization**: Internal graph analysis to minimize contraction cost.

## Performance Characteristics
- **Speed**: Comparable to ITensor in benchmarks; overhead is minimal.
- **Optimization**: "Behind-the-scenes" memory management.

## Comparison with Other Codes
- **vs ITensor**: Very similar philosophy (intelligent indices, C++). ITensor (Julia) is now more popular; Uni10 remains a strong C++ option with Python bindings.
- **vs TeNPy**: TeNPy is pure Python (mostly) + C modules; Uni10 is pure C++ + Python bindings.
- **vs TensorNetwork**: Uni10 is a specialized physics library; TensorNetwork is a general ML/Physics library.

## Application Areas
- **Spin Chains**: Heisenberg, Hubbard models (1D).
- **2D Systems**: PEPS simulations of square lattices.
- **Critical Systems**: MERA simulations of quantum critical points.

## Best Practices
- **Prototyping**: Use `pyUni10` to verify algorithms on small systems before moving to C++ for production.
- **Symmetries**: Always use `Uni10::UniTensor` with symmetries enabled for significant speedups.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/yingjerkao/uni10
2. Paper: "Uni10: An open-source library for tensor network algorithms" (arXiv:1511.05436).

**Confidence**: VERIFIED - Published and maintained.
**Verification status**: ✅ VERIFIED
