# DQC

## Official Resources
- Homepage: https://github.com/diffqc/dqc
- Documentation: https://github.com/diffqc/dqc#readme
- Source Repository: https://github.com/diffqc/dqc
- License: Apache License 2.0

## Overview
DQC (Differentiable Quantum Chemistry) is an open-source Python simulation code using PyTorch and xitorch for differentiable DFT and Hartree-Fock calculations. It enables automatic differentiation through quantum chemistry calculations, facilitating machine learning applications and gradient-based optimization.

**Scientific domain**: Differentiable quantum chemistry, ML/QC integration  
**Target user community**: Researchers combining quantum chemistry with machine learning

## Theoretical Methods
- Hartree-Fock (RHF, UHF)
- Density Functional Theory
- LDA, GGA functionals
- Orbital-free DFT
- Self-consistent field methods
- Differentiable implementations

## Capabilities (CRITICAL)
- PyTorch-based implementation
- Automatic differentiation through SCF
- GPU acceleration via PyTorch
- End-to-end gradients
- Neural network integration
- Learnable functionals
- Differentiable forces
- Batched calculations
- Custom loss functions
- Gradient-based optimization

## Key Strengths

### Differentiability:
- Full automatic differentiation
- Gradients through SCF
- Backpropagation support
- Custom objectives
- End-to-end training

### PyTorch Integration:
- Native tensor operations
- GPU support
- Neural network layers
- Optimizer compatibility
- Ecosystem integration

### ML Applications:
- Neural network potentials
- Learnable XC functionals
- Property prediction
- Inverse design
- Active learning

### Flexibility:
- Custom Hamiltonians
- Orbital-free DFT
- Novel approximations
- Research platform

## Inputs & Outputs
- **Input formats**:
  - PyTorch tensors
  - Molecular specifications
  - Python API
  
- **Output data types**:
  - Differentiable energies
  - Forces (automatic)
  - Densities
  - Gradients

## Interfaces & Ecosystem
- **PyTorch**: Native integration
- **xitorch**: Differentiable linear algebra
- **NumPy**: Array compatibility
- **ML frameworks**: Training pipelines

## Advanced Features

### Differentiable SCF:
- Implicit differentiation
- Stable gradients
- Convergence handling
- Adjoint methods

### Neural Functionals:
- Learnable exchange-correlation
- Neural network architectures
- Physics constraints
- Training procedures

### Optimization:
- Geometry optimization
- Property optimization
- Constrained optimization
- Multi-objective

## Performance Characteristics
- **Speed**: PyTorch GPU acceleration
- **Accuracy**: Standard DFT accuracy
- **System size**: Small-medium molecules
- **Memory**: PyTorch memory model
- **Parallelization**: GPU via PyTorch

## Computational Cost
- **HF/DFT**: Efficient on GPU
- **Differentiation**: Moderate overhead
- **Training**: Depends on model
- **Batching**: Efficient
- **Typical**: Suitable for ML workflows

## Limitations & Known Constraints
- **Method scope**: HF/DFT focus
- **System size**: Best for smaller molecules
- **Production**: More research-oriented
- **Documentation**: Basic
- **Community**: Small but active
- **Features**: Limited vs production codes

## Comparison with Other Codes
- **vs MESS**: Both differentiable; DQC PyTorch, MESS JAX
- **vs PySCF**: DQC differentiable-first, PySCF general
- **vs GPU4PySCF**: DQC PyTorch, GPU4PySCF CUDA
- **vs TorchANI**: DQC quantum chemistry, TorchANI potentials
- **Unique strength**: PyTorch-native differentiable QC

## Application Areas

### Machine Learning:
- Neural network training
- Differentiable simulations
- Property prediction
- Molecular design

### Functional Development:
- Learnable XC functionals
- Hybrid approaches
- Data-driven methods
- Novel approximations

### Optimization:
- Geometry optimization via gradients
- Property optimization
- Inverse problems
- Constrained design

## Best Practices

### PyTorch Usage:
- Understand autograd
- Efficient tensor operations
- GPU memory monitoring
- Gradient accumulation

### Training:
- Appropriate architectures
- Regularization
- Validation strategies
- Physical constraints

## Community and Support
- Open-source Apache 2.0
- GitHub development
- Related publications
- Academic collaborations
- Growing community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/diffqc/dqc
2. xitorch library: https://github.com/xitorch/xitorch
3. Academic publications on differentiable QC
4. PyTorch ecosystem

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, Apache 2.0)
- Documentation: README and examples
- Active development: Yes
- Research applications: Growing
