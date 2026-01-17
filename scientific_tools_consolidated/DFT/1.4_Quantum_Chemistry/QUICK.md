# QUICK

## Official Resources
- Homepage: https://quick-docs.readthedocs.io/
- Documentation: https://quick-docs.readthedocs.io/en/latest/
- Source Repository: https://github.com/merzlab/QUICK
- License: Mozilla Public License 2.0

## Overview
QUICK (QUantum Interaction Computational Kernel) is a GPU-enabled ab initio and density functional theory software package developed by the Götz and Merz labs. It leverages CUDA for GPU acceleration, providing significant speedups for electronic structure calculations on modern HPC hardware.

**Scientific domain**: Molecular quantum chemistry, GPU-accelerated calculations  
**Target user community**: Researchers needing fast HF/DFT calculations for medium-to-large molecules on GPU clusters

## Theoretical Methods
- Hartree-Fock (RHF, UHF)
- Density Functional Theory (pure and hybrid)
- LDA, GGA, meta-GGA, hybrid functionals
- Dispersion corrections (DFT-D3)
- COSMO implicit solvation
- Geometry optimization
- Molecular dynamics
- QM/MM (Amber interface)

## Capabilities (CRITICAL)
- GPU-accelerated two-electron integrals
- CUDA implementation of ERI evaluation
- Hybrid CPU/GPU execution
- Large-scale molecular calculations
- Direct SCF algorithm
- Cutoffs for linear scaling
- Geometry optimization
- Amber/QUICK QM/MM interface
- MPI multi-GPU parallelism
- Single and double precision

## Key Strengths

### GPU Acceleration:
- CUDA-optimized ERI kernels
- 10-100x speedups over CPU
- Multi-GPU scaling
- Modern NVIDIA GPU support
- Mixed precision capabilities

### HF/DFT Performance:
- Efficient integral evaluation
- Schwarz screening
- Grid-based DFT
- Hybrid functional support
- Large molecule capability

### Amber Integration:
- Native QM/MM interface
- Biochemical applications
- Enzyme active sites
- Drug-receptor interactions

### Scalability:
- Multi-GPU support
- MPI parallelism
- Linear-scaling techniques
- HPC cluster deployment

## Inputs & Outputs
- **Input formats**:
  - QUICK input files
  - Coordinate specifications
  - Basis set files
  - Amber interface files
  
- **Output data types**:
  - Energies and gradients
  - Molecular orbitals
  - Population analysis
  - Optimized structures
  - QM/MM energies

## Interfaces & Ecosystem
- **MD integration**: Amber molecular dynamics
- **Basis sets**: Standard Gaussian basis sets
- **Libraries**: CUDA, cuBLAS, MPI
- **Visualization**: Standard formats

## Advanced Features

### GPU Optimization:
- Custom CUDA kernels
- Batched integral evaluation
- Memory management
- Thread coarsening
- Shared memory usage

### QM/MM Calculations:
- Electrostatic embedding
- Mechanical embedding
- Link atom scheme
- Amber force fields

### Linear Scaling:
- Schwarz inequality screening
- Distance cutoffs
- Sparse matrix techniques
- Large system efficiency

## Performance Characteristics
- **Speed**: 10-100x GPU acceleration
- **Accuracy**: Standard HF/DFT accuracy
- **System size**: Hundreds to thousands of atoms
- **Memory**: GPU memory dependent
- **Parallelization**: Multi-GPU MPI scaling

## Computational Cost
- **HF**: Highly efficient on GPU
- **DFT**: Fast with grid optimization
- **Hybrid DFT**: GPU-accelerated exchange
- **Gradient**: Efficient for optimization
- **Typical**: Orders of magnitude faster than CPU

## Limitations & Known Constraints
- **GPU requirement**: NVIDIA CUDA GPUs needed
- **Method scope**: HF/DFT focus (no post-HF)
- **Basis sets**: Some limitations
- **Memory**: Limited by GPU memory
- **Features**: Fewer than general-purpose codes
- **Relativistic**: Not supported

## Comparison with Other Codes
- **vs TeraChem**: Both GPU-accelerated; different focus
- **vs Gaussian**: QUICK much faster on GPU, fewer methods
- **vs GAMESS**: QUICK GPU-native, GAMESS more features
- **vs PySCF-GPU4**: Different implementation approaches
- **Unique strength**: CUDA-optimized ERI, Amber QM/MM integration

## Application Areas

### Biochemistry:
- Enzyme mechanisms
- Drug binding energies
- Protein-ligand interactions
- Large biomolecular QM/MM

### Large Molecules:
- Organic semiconductors
- Polymers
- Supramolecular systems
- Nanoscale materials

### High-Throughput:
- Virtual screening
- Property prediction
- Database generation
- ML training data

## Best Practices

### GPU Utilization:
- Use modern NVIDIA GPUs
- Appropriate batch sizes
- Memory management
- Multi-GPU for large systems

### Calculation Setup:
- Appropriate cutoffs
- Grid density for DFT
- Convergence criteria
- Basis set selection

## Community and Support
- Open-source MPL 2.0
- Active GitHub development
- Merz and Götz research groups
- Documentation and tutorials
- Academic publications

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/merzlab/QUICK
2. Documentation: https://quick-docs.readthedocs.io/
3. Manathunga et al., J. Chem. Theory Comput. (2020) - QUICK paper
4. Active development at Michigan State and SDSC

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, MPL 2.0)
- Documentation: Comprehensive
- Active development: Yes
- Academic citations: Growing
