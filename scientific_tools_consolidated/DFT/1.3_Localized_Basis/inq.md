# inq

## Official Resources
- Homepage: https://gitlab.com/npneq/inq
- Documentation: https://inq.readthedocs.io/
- Source Repository: https://gitlab.com/npneq/inq (also GitHub mirror)
- GitHub: https://github.com/LLNL/inq
- License: Mozilla Public License 2.0

## Overview
inq is a modern, GPU-accelerated electronic structure code for DFT and real-time TDDFT calculations, developed at Lawrence Livermore National Laboratory (LLNL). Written in C++, it is designed from scratch for GPU execution and focuses on nonequilibrium phenomena in materials.

**Scientific domain**: Real-time dynamics, excited states, nonequilibrium phenomena, materials science  
**Target user community**: Researchers studying ultrafast dynamics, electronic excitations, and nonequilibrium materials

## Theoretical Methods
- Density Functional Theory (DFT)
- Real-time TDDFT
- Plane-wave/real-space hybrid basis
- LDA, GGA, meta-GGA functionals
- Hybrid functionals (HSE, PBE0)
- Spin-polarized and non-collinear spin
- Norm-conserving pseudopotentials

## Capabilities (CRITICAL)
- Ground-state DFT
- Real-time TDDFT simulations
- GPU-accelerated execution
- Absorption spectra
- Charge dynamics
- Ion dynamics
- Hybrid functionals
- Non-collinear magnetism
- Multiple spin configurations
- HPC cluster support

**Sources**: LLNL, arXiv publications, GitHub

## Key Strengths

### GPU Native:
- Designed for GPU from ground up
- CUDA/ROCm support
- Excellent GPU utilization
- Modern HPC architecture

### Real-Time TDDFT:
- Nonperturbative dynamics
- Ultrafast phenomena
- Charge carrier dynamics
- Strong-field physics

### Modern C++:
- Clean codebase (~12,000 lines)
- Maintainable
- Extensible
- Modern practices

### LLNL Quality:
- DOE national lab development
- Exascale computing focus
- Production-tested

## Inputs & Outputs
- **Input formats**:
  - C++ API
  - Python interface (developing)
  - Structure inputs
  
- **Output data types**:
  - Energies
  - Time-dependent observables
  - Spectra
  - Dynamics trajectories

## Interfaces & Ecosystem
- **HPC integration**:
  - MPI parallelization
  - GPU clusters
  - Slurm compatible
  
- **NPNEQ Center**:
  - Part of DOE center
  - Collaborative development

## Advanced Features

### Ultrafast Dynamics:
- Femtosecond timescales
- Pump-probe simulations
- Electron-ion coupling
- Hot carrier dynamics

### Strong-Field Physics:
- High-intensity excitations
- Nonlinear response
- Field-induced phenomena

### Non-Collinear Spin:
- Spin dynamics
- Magnetic systems
- Spin-orbit effects

## Performance Characteristics
- **Speed**: GPU-optimized
- **Accuracy**: Standard DFT/TDDFT
- **System size**: HPC-scalable
- **Memory**: GPU memory constraints
- **Parallelization**: MPI + GPU

## Computational Cost
- **GPU efficiency**: Excellent
- **RT-TDDFT**: Time-stepping overhead
- **Typical**: GPU cluster runs

## Limitations & Known Constraints
- **Maturity**: Actively developing
- **Documentation**: Growing
- **Traditional DFT**: TDDFT focus
- **CPU fallback**: GPU-primary

## Comparison with Other Codes
- **vs Octopus**: Both RT-TDDFT, inq GPU-native
- **vs SALMON**: Similar TDDFT, different implementation
- **Unique strength**: GPU-native design, LLNL backing

## Application Areas

### Ultrafast Science:
- Pump-probe spectroscopy
- Attosecond dynamics
- Carrier relaxation

### Optical Properties:
- Absorption spectra
- Dielectric response
- Plasmonics

### Materials Dynamics:
- Electron-phonon coupling
- Phase transitions
- Radiation effects
- Energy transfer

## Best Practices

### Ground State Setup:
- Converge ground state fully before dynamics
- Use appropriate exchange-correlation functional
- Test k-point convergence for periodic systems
- Verify energy conservation

### RT-TDDFT Simulations:
- Choose appropriate time step (0.01-0.1 fs typical)
- Apply kick pulse for linear response spectra
- Use proper absorbing boundaries for finite systems
- Monitor total energy drift

### GPU Optimization:
- Match problem size to GPU memory
- Use batched calculations when possible
- Profile to find memory bottlenecks
- Consider multi-GPU for larger systems

### Post-Processing:
- Fourier transform for spectra from time-domain data
- Apply windowing functions to reduce spectral artifacts
- Calculate dipole moments for absorption spectra
- Analyze induced densities for insight

## Computational Cost
- **Ground state DFT**: Comparable to other GPU codes
- **RT-TDDFT**: Linear in simulation time, each step similar cost to ground state
- **GPU efficiency**: 10-100x speedup over CPU for suitable systems
- **Memory**: GPU memory often limiting factor
- **Typical**: Seconds for small molecules, hours for large periodic
- **Scaling**: Cubic O(N³) with system size for ground state

## Community and Support
- Open source MPL 2.0
- LLNL development
- DOE NPNEQ center
- GitLab/GitHub
- Active development

## Verification & Sources
**Primary sources**:
1. GitLab: https://gitlab.com/npneq/inq
2. GitHub: https://github.com/LLNL/inq
3. arXiv publications
4. LLNL computational materials

**Confidence**: VERIFIED - LLNL official code

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MPL 2.0)
- Documentation: ReadTheDocs
- Development: Active
- Specialty: GPU-native DFT/RT-TDDFT, nonequilibrium dynamics
