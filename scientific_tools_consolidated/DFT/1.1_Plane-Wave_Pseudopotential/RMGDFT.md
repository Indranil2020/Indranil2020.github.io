# RMG (Real-space Multigrid DFT)

## Official Resources
- Homepage: https://github.com/RMGDFT/rmgdft
- Documentation: https://github.com/RMGDFT/rmgdft/blob/master/README.md
- Source Repository: https://github.com/RMGDFT/rmgdft
- License: GNU General Public License v2.0

## Overview
RMG (Real-space Multigrid) is an open-source DFT code using real-space multigrid methods for solving the Kohn-Sham equations. Developed at North Carolina State University, RMG employs multigrid preconditioning for efficient convergence and is designed for excellent parallel performance on modern HPC systems with GPU acceleration. It is particularly effective for large-scale electronic structure calculations and materials simulations requiring high accuracy.

**Scientific domain**: Real-space DFT, multigrid methods, GPU computing, large systems  
**Target user community**: Materials scientists, HPC users, large-system researchers

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA)
- Real-space multigrid representation
- Finite-difference discretization
- Norm-conserving pseudopotentials
- Ultrasoft pseudopotentials
- Multigrid preconditioning
- van der Waals corrections
- DFT+U for correlated systems
- Spin-orbit coupling
- Non-collinear magnetism

## Capabilities (CRITICAL)
- Ground state electronic structure
- Geometry optimization
- Molecular dynamics (NVE, NVT)
- Band structures and DOS
- Forces and stress tensors
- Large system calculations
- Multigrid acceleration
- GPU acceleration (CUDA)
- Excellent parallel scalability
- Real-space representation
- Periodic and non-periodic systems
- Efficient convergence
- High accuracy

**Sources**: GitHub repository (https://github.com/RealSpaceGroup/RMGDFT)

## Key Strengths

### Multigrid Methods:
- Fast convergence
- Efficient preconditioning
- Hierarchical grids
- Optimal scaling
- Reduced iterations

### GPU Acceleration:
- Extensive CUDA support
- Dramatic speedups
- Multi-GPU capable
- Optimized kernels
- Production-ready

### Real-Space:
- Natural for non-periodic
- No FFT overhead
- Local operations
- Sparse matrices
- Scalability benefits

### Scalability:
- Excellent parallel efficiency
- MPI+GPU hybrid
- Large-scale systems
- HPC optimized

### Open Source:
- GPL v2 licensed
- Active GitHub development
- Community contributions
- Well-maintained

## Inputs & Outputs
- **Input formats**:
  - Text-based input file
  - XYZ coordinates
  - Pseudopotential files
  - Cell parameters
  
- **Output data types**:
  - Text output
  - Wavefunctions
  - Charge densities
  - Energy and forces
  - Trajectory files

## Interfaces & Ecosystem
- **Preparation**:
  - Standard format conversion
  - Python scripts
  - ASE integration (developing)
  
- **Analysis**:
  - Custom tools
  - Python post-processing
  - Standard viewers
  
- **Pseudopotentials**:
  - Norm-conserving
  - Ultrasoft
  - Standard formats
  
- **Parallelization**:
  - MPI parallelization
  - GPU offloading (CUDA)
  - Hybrid MPI+GPU
  - Domain decomposition

## Workflow and Usage

### Example Input:

```
# Silicon calculation
start_mode = "LCAO Start"
calculation_mode = "Quench Electrons"

latticevec = "
  5.13 0.0 0.0
  0.0 5.13 0.0
  0.0 0.0 5.13
"

atoms = "
  Si 0.0 0.0 0.0 1 1 1
  Si 0.25 0.25 0.25 1 1 1
"

pseudopotential = "Si.UPF"

wavefunction_grid = "64 64 64"
kpoint_mesh = "4 4 4"

xc_type = "GGA PBE"
```

### Running RMG:
```bash
rmg input.in
mpirun -np 8 rmg input.in
# GPU
mpirun -np 8 rmg-gpu input.in
```

## Advanced Features

### Multigrid Preconditioning:
- Hierarchical grid levels
- V-cycle or W-cycle
- Fast convergence
- Reduced SCF iterations
- Optimal complexity

### GPU Acceleration:
- CUDA implementation
- Significant speedup (5-10x)
- Multi-GPU scaling
- Optimized for NVIDIA
- Memory efficient

### Real-Space Grids:
- Finite differences
- Direct diagonalization
- Sparse operations
- Natural boundaries
- Systematic convergence

### Parallel Performance:
- MPI domain decomposition
- GPU parallelization
- Hybrid approach
- Good scaling
- HPC ready

## Performance Characteristics
- **Speed**: Very fast with GPU
- **Scaling**: Good parallel scaling
- **GPU**: 5-10x acceleration
- **Memory**: Moderate requirements
- **Typical systems**: 100-1000 atoms

## Computational Cost
- **DFT**: Efficient
- **Multigrid**: Faster convergence
- **GPU**: Dramatically reduced time
- **Large systems**: Feasible
- **MD**: Production runs possible

## Limitations & Known Constraints
- **Smaller community**: Less established
- **Documentation**: Basic
- **Features**: Fewer than major codes
- **Learning curve**: Moderate
- **GPU**: CUDA only (NVIDIA)
- **Platform**: Linux primarily

## Comparison with Other Codes
- **vs VASP**: RMG open-source, real-space multigrid
- **vs Quantum ESPRESSO**: RMG multigrid acceleration, GPU
- **vs PARSEC**: Both real-space, RMG has multigrid
- **vs SPARC**: Similar real-space approach
- **Unique strength**: Multigrid methods, GPU acceleration, open-source

## Application Areas

### Large Systems:
- Nanostructures
- Complex materials
- Interfaces
- Amorphous systems

### Materials Science:
- Electronic structure
- Structural properties
- Phase stability
- Defects

### GPU Computing:
- HPC applications
- Fast turnaround
- Large-scale screening
- Production calculations

## Best Practices

### Grid Convergence:
- Test wavefunction grid
- Check charge grid
- Systematic refinement
- Balance accuracy/cost

### Multigrid Setup:
- Appropriate grid levels
- V-cycle standard
- Monitor convergence
- Optimize cycles

### GPU Usage:
- Use GPU build
- Balance CPU-GPU
- Optimize batch sizes
- Monitor utilization

### Parallelization:
- Test scaling
- Optimize MPI layout
- Balance domains
- Use GPU when available

## Community and Support
- Open-source (GPL v2)
- GitHub repository
- Issue tracking
- Basic documentation
- Academic development

## Educational Resources
- README on GitHub
- Example inputs
- Published papers
- Source code documentation

## Development
- Active GitHub
- NC State University
- Community contributions
- Regular updates
- Modern codebase

## Research Applications
- Large-scale DFT
- Materials discovery
- Electronic structure
- Method development

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/RealSpaceGroup/RMGDFT
2. README: https://github.com/RealSpaceGroup/RMGDFT/blob/master/README.md
3. E. L. Briggs et al., Phys. Rev. B 54, 14362 (1996) - Multigrid methods for DFT

**Secondary sources**:
1. GitHub documentation
2. Published studies using RMG
3. Real-space multigrid literature
4. HPC conference proceedings

**Confidence**: VERIFIED - GitHub repository confirmed

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: Basic (README, source)
- Source code: OPEN (GitHub, GPL v2)
- Community support: GitHub issues
- Active development: Regular commits
- Specialized strength: Real-space multigrid methods, GPU acceleration, efficient convergence, open-source
