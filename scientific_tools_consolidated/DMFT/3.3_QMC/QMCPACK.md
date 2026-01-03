# QMCPACK (Quantum Monte Carlo Package)

## Official Resources
- Homepage: https://qmcpack.org/
- Documentation: https://qmcpack.readthedocs.io/
- Source Repository: https://github.com/QMCPACK/qmcpack
- License: BSD 3-Clause License (open-source)

## Overview
QMCPACK is a modern, high-performance implementation of continuum quantum Monte Carlo (QMC) methods for electronic structure calculations of molecules, 2D materials, and solids. Developed as a community code with major contributions from Oak Ridge National Laboratory, QMCPACK implements Variational Monte Carlo (VMC), Diffusion Monte Carlo (DMC), and related methods. It is optimized for leadership-class supercomputers and provides production-quality calculations for realistic materials with unprecedented accuracy.

**Scientific domain**: Quantum Monte Carlo, electronic structure, materials science  
**Target user community**: Materials scientists, computational chemists, HPC users

## Theoretical Methods
- Variational Monte Carlo (VMC)
- Diffusion Monte Carlo (DMC)
- Reptation QMC (RQMC)
- Fixed-node approximation
- Slater-Jastrow wavefunctions
- Multi-determinant expansions
- Pseudopotentials
- All-electron calculations
- Real-space methods

## Capabilities (CRITICAL)
**Category**: Open-source QMC code
- VMC and DMC methods
- Molecules and solids
- Periodic boundary conditions
- Slater-Jastrow wavefunctions
- Multi-reference trial functions
- Backflow correlations
- GPU acceleration (CUDA, HIP)
- MPI + OpenMP + GPU
- Wavefunction optimization
- Energy calculations
- Forces and structural optimization
- Excited states
- DFT trial wavefunction input
- AFQMC (Auxiliary-Field QMC)
- Production quality

**Sources**: Official website, documentation, publications

## Key Strengths

### Accuracy:
- Beyond-DFT precision
- Benchmark-quality results
- Systematic improvement
- Fixed-node DMC
- Chemical accuracy achievable

### HPC Performance:
- Leadership-class scaling
- GPU acceleration
- Hybrid parallelization
- Optimized kernels
- Exascale-ready

### Production Quality:
- Well-tested
- Large community
- Active development
- Comprehensive documentation
- Publication-ready

### Versatility:
- Molecules to solids
- Finite/periodic systems
- Various materials
- Multiple QMC methods
- Flexible workflows

## Inputs & Outputs
- **Input formats**:
  - XML input files
  - DFT trial wavefunctions (VASP, Quantum ESPRESSO, PySCF)
  - Pseudopotentials
  - Structure files
  
- **Output data types**:
  - Total energies
  - Forces
  - Structural properties
  - Excited states
  - Statistical data
  - HDF5 archives

## Interfaces & Ecosystem

### DFT Codes:
- Quantum ESPRESSO (pw2qmcpack)
- VASP
- PySCF
- Gaussian
- GAMESS

### Tools:
- Nexus workflow tool
- Python utilities
- Analysis scripts
- Wavefunction converters

## Workflow and Usage

### Installation:
```bash
# Clone repository
git clone https://github.com/QMCPACK/qmcpack.git
cd qmcpack
mkdir build && cd build

# Configure with GPU
cmake -DQMC_CUDA=1 ..
make -j8
```

### DFT Trial Wavefunction:
```bash
# Quantum ESPRESSO
pw.x < scf.in > scf.out
pw2qmcpack.x < p2q.in > p2q.out
```

### VMC Optimization:
```xml
<!-- opt.xml -->
<simulation>
  <qmc method="linear" move="pbyp">
    <parameter name="warmupSteps">100</parameter>
    <parameter name="blocks">200</parameter>
    <parameter name="timestep">0.5</parameter>
    <parameter name="samples">16000</parameter>
  </qmc>
</simulation>
```

### DMC Calculation:
```xml
<!-- dmc.xml -->
<simulation>
  <qmc method="dmc" move="pbyp">
    <parameter name="targetWalkers">1920</parameter>
    <parameter name="blocks">1000</parameter>
    <parameter name="timestep">0.01</parameter>
    <parameter name="warmupSteps">50</parameter>
  </qmc>
</simulation>
```

### Run QMCPACK:
```bash
# MPI + GPU
mpirun -n 8 qmcpack opt.xml
mpirun -n 8 qmcpack dmc.xml
```

### Using Nexus:
```python
from nexus import settings, generate_physical_system
from nexus import generate_pwscf, generate_pw2qmcpack
from nexus import generate_qmcpack, loop, run_project

# Setup system
system = generate_physical_system(
    structure='diamond.xsf',
    C=4
)

# DFT calculation
scf = generate_pwscf(
    system=system,
    job=job_scf,
    input_dft='lda'
)

# Convert to QMCPACK
p2q = generate_pw2qmcpack(
    dependencies=(scf,'orbitals')
)

# VMC optimization
opt = generate_qmcpack(
    system=system,
    job=job_opt,
    input_type='basic',
    qmc='opt',
    dependencies=(p2q,'orbitals')
)

# DMC
dmc = generate_qmcpack(
    system=system,
    job=job_dmc,
    input_type='basic',
    qmc='dmc',
    dependencies=(opt,'jastrow')
)

run_project()
```

## Advanced Features

### GPU Acceleration:
- CUDA support
- HIP (AMD ROCm)
- Mixed precision
- Performance portable
- Exascale systems

### AFQMC:
- Auxiliary-Field QMC
- Finite temperature
- Correlated systems
- Complementary to VMC/DMC

### Wavefunction Optimization:
- Linear method
- Energy minimization
- Variance minimization
- Efficient optimization

### Excited States:
- Excited state calculations
- Promotion methods
- Optical gaps
- Multiple states

## Performance Characteristics
- **Speed**: HPC-optimized, GPU-accelerated
- **Accuracy**: Chemical accuracy achievable
- **Scalability**: Exascale-ready
- **System size**: 100s-1000s electrons
- **Typical**: Leadership-class HPC

## Computational Cost
- Expensive but accurate
- DMC: O(N³-N⁴) electrons
- GPU acceleration crucial
- HPC resources required
- Production: Millions of core-hours

## Limitations & Known Constraints
- **Fixed-node approximation**: DMC nodal error
- **Computational cost**: Very expensive
- **HPC required**: Not for desktop
- **Trial wavefunction**: Quality matters
- **Learning curve**: Steep
- **Pseudopotentials**: Careful selection needed

## Comparison with Other QMC Codes
- **vs CASINO**: QMCPACK HPC-focused, CASINO feature-rich
- **vs TurboRVB**: QMCPACK larger community, TurboRVB specialized
- **Unique strength**: HPC performance, GPU acceleration, community support, production quality, exascale-ready

## Application Areas

### Materials Science:
- Crystals and solids
- 2D materials
- Defects
- Surfaces
- Nanostructures

### Chemistry:
- Molecular systems
- Reaction energies
- Excited states
- Chemical accuracy

### Condensed Matter:
- Electronic structure
- Benchmark calculations
- Beyond-DFT accuracy
- Correlation effects

### High-Performance Computing:
- Exascale applications
- GPU computing
- Performance optimization
- Scalability studies

## Best Practices

### Trial Wavefunctions:
- Quality DFT starting point
- Hybrid functionals preferred
- Multi-determinant when needed
- Jastrow optimization crucial

### DMC Parameters:
- Timestep extrapolation
- Population control
- Finite-size corrections
- Error analysis

### HPC Usage:
- GPU acceleration
- Load balancing
- I/O optimization
- Checkpoint/restart

## Community and Support
- Open-source (BSD 3-Clause)
- Large user community
- Active development
- Annual workshops
- Mailing lists
- Slack channel
- GitHub issue tracking

## Educational Resources
- Comprehensive documentation
- Tutorials
- Workshop materials
- Example inputs
- Publication list
- QMC schools

## Development
- ORNL leadership
- Multi-institutional
- Community contributions
- Active development
- Exascale Computing Project
- Regular releases

## Research Impact
QMCPACK is used for high-precision electronic structure calculations across chemistry, materials science, and condensed matter physics, with thousands of publications and major allocations on leadership-class supercomputers worldwide.

## Verification & Sources
**Primary sources**:
1. Homepage: https://qmcpack.org/
2. Documentation: https://qmcpack.readthedocs.io/
3. GitHub: https://github.com/QMCPACK/qmcpack
4. Publications: J. Chem. Theory Comput. 16, 4 (2020)

**Secondary sources**:
1. User publications
2. QMC literature
3. HPC reports

**Confidence**: CONFIRMED - Leading QMC code

**Verification status**: ✅ CONFIRMED
- Website: ACTIVE
- License: BSD 3-Clause (open-source)
- **Category**: Open-source QMC code
- Status: Actively developed
- Community: Very large, international
- Specialized strength: High-performance quantum Monte Carlo for electronic structure, GPU acceleration, exascale-ready, VMC/DMC/AFQMC methods, production quality, leadership-class HPC, benchmark accuracy, large community, comprehensive documentation, materials and molecules
