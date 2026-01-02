# ABACUS

## Official Resources
- Homepage: https://abacus.deepmodeling.com/
- Documentation: https://abacus.deepmodeling.com/en/latest/
- Source Repository: https://github.com/deepmodeling/abacus-develop
- License: GNU Lesser General Public License v3.0

## Overview
ABACUS (Atomic-orbital Based Ab-initio Computation at UStc) is an open-source DFT package supporting both plane-wave and numerical atomic orbital (NAO) basis sets. Developed by the DeepModeling community in China, ABACUS is designed for efficient electronic structure calculations, materials simulations, and integration with machine learning workflows. It features excellent GPU acceleration, hybrid basis sets, and is particularly strong in Chinese academic and industrial applications.

**Scientific domain**: DFT with multiple basis sets, materials science, machine learning integration  
**Target user community**: Materials scientists, ML researchers, Chinese research community

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Plane-wave basis sets
- Numerical atomic orbitals (NAO)
- Hybrid plane-wave/NAO basis
- Norm-conserving pseudopotentials
- Stochastic DFT methods
- van der Waals corrections (DFT-D2/D3)
- DFT+U for correlated systems
- Hybrid functionals
- Spin-orbit coupling
- Non-collinear magnetism

## Capabilities (CRITICAL)
- Ground state electronic structure
- Geometry optimization
- Molecular dynamics (Born-Oppenheimer)
- Cell relaxation
- Band structures and DOS
- Phonons via DFPT
- Elastic constants
- Surface calculations
- Machine learning interface
- DeepKS (ML correction)
- GPU acceleration
- Multiple basis sets (PW, NAO, hybrid)
- Stochastic DFT for large systems
- Linear-scaling algorithms
- Excellent parallelization
- Python interface (PyABACUS)

**Sources**: Official ABACUS documentation (https://abacus.deepmodeling.com/), GitHub repository

## Key Strengths

### Multiple Basis Sets:
- Plane-waves (PW)
- Numerical atomic orbitals (NAO)
- Hybrid PW/NAO
- Flexibility in choice
- Method comparison

### Machine Learning:
- DeepModeling integration
- DeepKS correction
- ML potential interface
- DP-GEN workflows
- Modern AI integration

### GPU Acceleration:
- CUDA support
- Significant speedup
- Production-ready
- Optimized kernels

### Stochastic DFT:
- Large system capability
- Linear scaling
- Statistical sampling
- Thousands of atoms

### Open Source:
- LGPL v3 licensed
- Active GitHub development
- Growing community
- Chinese and international

## Inputs & Outputs
- **Input formats**:
  - INPUT file (main parameters)
  - STRU file (structure)
  - KPT file (k-points)
  - Pseudopotential files
  
- **Output data types**:
  - Text output (OUT.*)
  - Structure files
  - Wavefunction data
  - Band structure files
  - DOS files

## Interfaces & Ecosystem
- **DeepModeling**:
  - DP-GEN workflows
  - DeePMD integration
  - DeepKS correction
  - ML ecosystem
  
- **Analysis**:
  - PyABACUS (Python interface)
  - Custom tools
  - Standard visualization
  
- **Parallelization**:
  - MPI parallelization
  - GPU acceleration (CUDA)
  - Hybrid parallelization
  - Good scaling

## Workflow and Usage

### Example INPUT File:

```
INPUT_PARAMETERS
calculation  scf
basis_type   lcao

ecutwfc      100
scf_thr      1.0e-7
scf_nmax     100

smearing_method  gaussian
smearing_sigma   0.001

kspacing     0.1
```

### Example STRU File:

```
ATOMIC_SPECIES
Si 28.0855 Si_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
5.43

LATTICE_VECTORS
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5

ATOMIC_POSITIONS
Direct

Si
0.0
2
0.00 0.00 0.00 1 1 1
0.25 0.25 0.25 1 1 1
```

### Running ABACUS:
```bash
abacus
mpirun -np 16 abacus
```

## Advanced Features

### Stochastic DFT:
- Statistical sampling
- Large systems (1000+ atoms)
- Linear-scaling approach
- Reduced computational cost
- Controlled accuracy

### DeepKS:
- Machine learning correction
- Improved accuracy
- Trained models
- DeepModeling framework
- Enhanced DFT

### Multiple Basis Options:
- Choose PW for accuracy
- Choose NAO for efficiency
- Hybrid for balance
- System-dependent optimization

### GPU Acceleration:
- Offload to GPU
- CUDA kernels
- Speedup factors
- Multi-GPU support

### DFPT:
- Phonon calculations
- Dielectric response
- Born effective charges
- Efficient implementation

## Performance Characteristics
- **Speed**: Competitive
- **Scaling**: Good parallel performance
- **GPU**: Significant acceleration
- **Memory**: Moderate
- **Typical systems**: 50-1000 atoms

## Computational Cost
- **DFT (PW)**: Standard
- **DFT (NAO)**: More efficient
- **Stochastic**: Very efficient for large
- **GPU**: Dramatically faster
- **ML corrections**: Minimal overhead

## Limitations & Known Constraints
- **Smaller community**: Less established globally
- **Documentation**: Good but primarily Chinese/English
- **Learning curve**: Moderate
- **Features**: Fewer than VASP/QE
- **Platform**: Linux primarily
- **GPU**: CUDA only

## Comparison with Other Codes
- **vs VASP**: ABACUS open-source, multiple basis
- **vs Quantum ESPRESSO**: Similar PW capabilities, ABACUS adds NAO
- **vs SIESTA**: Both NAO, ABACUS also PW
- **vs CP2K**: Different implementations
- **Unique strength**: Multiple basis sets, ML integration, stochastic DFT, DeepModeling ecosystem

## Application Areas

### Materials Science:
- Electronic structure
- Structural properties
- Phase stability
- Defects and dopants

### Machine Learning:
- Training data generation
- DeepKS applications
- ML potential development
- DP-GEN workflows

### Large Systems:
- Stochastic DFT
- Amorphous materials
- Complex systems
- Interfaces

### Method Development:
- Basis set comparison
- Algorithm testing
- ML integration

## Best Practices

### Basis Set Selection:
- PW for high accuracy
- NAO for large systems
- Test both for system
- Convergence testing

### Convergence:
- Plane-wave cutoff (PW)
- NAO basis size
- K-point sampling
- SCF threshold

### GPU Usage:
- Enable GPU acceleration
- Monitor utilization
- Balance CPU-GPU
- Optimize workload

### ML Integration:
- Use DeepKS when appropriate
- Train models carefully
- Validate corrections
- Follow DeepModeling best practices

## Community and Support
- Open-source (LGPL v3)
- Active GitHub repository
- Documentation website
- DeepModeling community
- Chinese and international users
- Regular updates

## Educational Resources
- Online documentation
- Tutorial examples
- GitHub wiki
- DeepModeling resources
- Published papers

## Development
- Active GitHub development
- DeepModeling project
- University of Science and Technology of China
- Community contributions
- Regular releases

## DeepModeling Ecosystem
- Part of DeepModeling community
- Integration with DeePMD
- DP-GEN workflows
- ML potential development
- AI for science

## Verification & Sources
**Primary sources**:
1. Official website: https://abacus.deepmodeling.com/
2. Documentation: https://abacus.deepmodeling.com/en/latest/
3. GitHub repository: https://github.com/deepmodeling/abacus-develop
4. M. Chen et al., Comput. Phys. Commun. 285, 108634 (2023) - ABACUS overview

**Secondary sources**:
1. ABACUS documentation and tutorials
2. DeepModeling community resources
3. Published studies using ABACUS
4. GitHub discussions

**Confidence**: VERIFIED - GitHub repository and official website confirmed

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, LGPL v3)
- Community support: GitHub issues, documentation
- Active development: Very active GitHub
- Specialized strength: Multiple basis sets (PW/NAO), ML integration, stochastic DFT, DeepModeling ecosystem, GPU acceleration
