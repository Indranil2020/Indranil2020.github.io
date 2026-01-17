# GPU4PySCF

## Official Resources
- Homepage: https://github.com/pyscf/gpu4pyscf
- Documentation: https://gpu4pyscf.readthedocs.io/
- Source Repository: https://github.com/pyscf/gpu4pyscf
- Parent Project: https://pyscf.org/
- License: Apache License 2.0

## Overview
GPU4PySCF is a GPU-accelerated extension for PySCF, providing CUDA implementations of two-electron repulsion integrals and DFT calculations. It enables significant speedups for Hartree-Fock and DFT calculations while maintaining compatibility with the PySCF ecosystem.

**Scientific domain**: GPU-accelerated quantum chemistry, HF/DFT  
**Target user community**: PySCF users needing GPU acceleration for larger molecules

## Theoretical Methods
- Restricted/Unrestricted Hartree-Fock
- Density Functional Theory
- LDA, GGA, meta-GGA, hybrid functionals
- Density fitting (DF/RI-J/K)
- Direct SCF
- Geometry optimization
- Hessian calculations

## Capabilities (CRITICAL)
- GPU-accelerated ERI evaluation
- CUDA kernel implementations
- Density fitting on GPU
- SCF acceleration
- DFT grid calculations
- Gradient calculations
- PySCF API compatibility
- Multi-GPU support
- cuBLAS integration
- Memory-efficient algorithms

## Key Strengths

### GPU Acceleration:
- CUDA-optimized integrals
- 10-100x speedups
- Modern GPU support
- Multi-GPU capability
- Memory management

### PySCF Integration:
- Drop-in replacement
- Same API
- Ecosystem compatibility
- Workflow preservation
- Easy adoption

### DFT Performance:
- Fast grid integration
- Hybrid functionals
- Exchange matrices
- Large molecules

### Density Fitting:
- RI-J and RI-K
- GPU-accelerated DF
- Auxiliary basis sets
- Efficient memory use

## Inputs & Outputs
- **Input formats**:
  - PySCF mol objects
  - Standard coordinates
  - Basis set strings
  
- **Output data types**:
  - Energies
  - Gradients
  - Densities
  - Orbitals
  - Properties

## Interfaces & Ecosystem
- **PySCF**: Full integration
- **NumPy/CuPy**: Array backends
- **CUDA**: GPU runtime
- **Post-HF**: Enable GPU-accelerated reference

## Advanced Features

### ERI Acceleration:
- Schwarz screening
- Shell pair sorting
- Memory blocking
- Batched evaluation

### Density Fitting GPU:
- Three-center integrals
- Coulomb/exchange
- Auxiliary screening
- Efficient contraction

### Gradient Calculations:
- Analytical gradients
- Geometry optimization
- Nuclear forces
- GPU-accelerated

## Performance Characteristics
- **Speed**: 10-100x GPU acceleration
- **Accuracy**: Identical to CPU PySCF
- **System size**: Significantly larger molecules
- **Memory**: GPU memory dependent
- **Parallelization**: Multi-GPU MPI

## Computational Cost
- **HF**: Dramatic GPU speedup
- **DFT**: Fast grid + integrals
- **Hybrid DFT**: Excellent acceleration
- **DF methods**: Very efficient
- **Typical**: Enables previously infeasible calculations

## Limitations & Known Constraints
- **GPU requirement**: NVIDIA CUDA needed
- **Memory**: Limited by GPU RAM
- **Methods**: HF/DFT focus (post-HF separate)
- **Basis sets**: Standard GTOs
- **Platform**: Linux primarily
- **Installation**: CUDA setup required

## Comparison with Other Codes
- **vs QUICK**: Both GPU; GPU4PySCF PySCF ecosystem
- **vs TeraChem**: GPU4PySCF open-source, TeraChem commercial
- **vs Psi4 GPU**: Different implementations
- **vs CPU PySCF**: Same results, much faster
- **Unique strength**: PySCF integration, open-source, community

## Application Areas

### Large Molecules:
- Proteins and enzymes
- Organic semiconductors
- Supramolecular systems
- Drug molecules

### High-Throughput:
- Database generation
- Virtual screening
- ML training data
- Property prediction

### Method Development:
- Reference calculations
- Benchmarking
- Algorithm testing
- Rapid prototyping

## Best Practices

### GPU Setup:
- Modern NVIDIA GPU
- Sufficient GPU memory
- CUDA toolkit installed
- cuBLAS optimized

### Calculation Strategy:
- Use density fitting
- Appropriate cutoffs
- Batch similar molecules
- Monitor GPU memory

## Community and Support
- Open-source Apache 2.0
- PySCF community support
- Active GitHub development
- Documentation and examples
- Academic publications

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/pyscf/gpu4pyscf
2. PySCF ecosystem: https://pyscf.org/
3. J. Chem. Theory Comput. (2024) GPU4PySCF paper
4. Active development

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, Apache 2.0)
- Documentation: ReadTheDocs
- Active development: Yes
- Part of PySCF project
