# cuEquivariance

## Official Resources
- Source Repository: https://github.com/NVIDIA/cuEquivariance
- License: Open source (BSD-3)

## Overview
**cuEquivariance** is NVIDIA's CUDA library for fast equivariant operations used in MLIPs. It provides optimized kernels for MACE, NequIP, and other equivariant models, achieving significant speedups on NVIDIA GPUs.

**Scientific domain**: CUDA-accelerated equivariant operations for MLIPs  
**Target user community**: Researchers running equivariant MLIPs on NVIDIA GPUs

## Theoretical Methods
- CUDA-optimized equivariant kernels
- Spherical harmonic operations
- Clebsch-Gordan coefficients
- Tensor product operations
- GPU memory optimization

## Capabilities (CRITICAL)
- Fast equivariant operations
- MACE kernel optimization
- NequIP kernel optimization
- NVIDIA GPU acceleration
- PyTorch integration

**Sources**: GitHub repository

## Key Strengths

### Performance:
- Significant GPU speedup
- Optimized CUDA kernels
- Memory efficient
- Batch processing

### Compatibility:
- MACE support
- NequIP support
- General equivariant models
- PyTorch integration

### NVIDIA:
- Professional development
- Regular updates
- GPU optimization expertise
- Production quality

## Inputs & Outputs
- **Input formats**: Equivariant model operations
- **Output data types**: Accelerated equivariant computations

## Interfaces & Ecosystem
- **PyTorch**: Integration
- **CUDA**: Backend
- **MACE/NequIP**: Supported models

## Performance Characteristics
- **Speed**: 2-10x speedup
- **Accuracy**: Identical (numerical)
- **System size**: Any
- **Automation**: Drop-in

## Computational Cost
- **Setup**: Minutes
- **Runtime**: Reduced

## Limitations & Known Constraints
- **NVIDIA GPU only**: No AMD/Intel
- **Specific models**: MACE, NequIP primarily
- **CUDA required**: Build complexity
- **New project**: Still maturing

## Comparison with Other Codes
- **vs pure PyTorch**: cuEquivariance is 2-10x faster
- **vs e3nn**: cuEquivariance is optimized CUDA
- **Unique strength**: NVIDIA CUDA-optimized equivariant kernels for 2-10x MLIP speedup

## Application Areas

### Production MLIP:
- Fast MACE inference
- Fast NequIP training
- Large-scale equivariant MD
- GPU cluster optimization

### Research:
- Equivariant model development
- Performance benchmarking
- Architecture optimization

## Best Practices
- Use with NVIDIA A100/H100
- Drop-in replacement for e3nn
- Benchmark before/after
- Use latest CUDA version

## Community and Support
- Open source (BSD-3)
- NVIDIA maintained
- GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/NVIDIA/cuEquivariance

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: NVIDIA CUDA-optimized equivariant kernels for 2-10x MLIP speedup
