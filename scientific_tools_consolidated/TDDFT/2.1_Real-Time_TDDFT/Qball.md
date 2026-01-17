# Qb@ll (Qball)

## Official Resources
- Homepage: https://computing.llnl.gov/projects/qball
- Documentation: https://github.com/LLNL/qball/blob/master/README.md
- Source Repository: https://github.com/LLNL/qball
- License: GNU General Public License v3.0 (LLNL-CODE-635376)

## Overview
Qb@ll (also written as Qball or qb@ll) is a first-principles molecular dynamics code developed at Lawrence Livermore National Laboratory. It computes electronic structure of atoms, molecules, solids, and liquids using Density Functional Theory with a plane-wave basis. Qb@ll is a fork of the Qbox code by François Gygi, optimized for high-performance computing including Real-Time TDDFT capabilities.

**Scientific domain**: Ab initio molecular dynamics, electronic structure, real-time electron dynamics  
**Target user community**: HPC users at national labs and institutions needing scalable plane-wave DFT/TDDFT

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis set
- Norm-conserving pseudopotentials
- DFT-GGA and Hybrid DFT functionals
- Born-Oppenheimer molecular dynamics
- Car-Parrinello molecular dynamics
- Real-Time TDDFT (RT-TDDFT) - developed in Qb@ch branch
- NVT simulations with stochastic thermostats

## Capabilities
- **Ground-state DFT** for all system types
- **First-principles molecular dynamics** (FPMD)
- **Born-Oppenheimer MD** with forces from DFT
- **Car-Parrinello MD** for efficient dynamics
- **Real-Time TDDFT** (via Qb@ch development)
- **Hybrid functionals** support
- **Large-scale parallel** execution
- **HPC optimized** (designed for supercomputers)

## Key Strengths

### HPC Performance:
- Designed for leadership-class supercomputers
- Excellent parallel scaling
- Blue Gene/Q optimized configurations
- MPI + OpenMP hybrid parallelization

### Qbox Foundation:
- Fork of established Qbox code
- Plane-wave accuracy
- Proven algorithms
- Active LLNL support

### RT-TDDFT Development:
- Qb@ch branch (UNC Chapel Hill) focuses on RT-TDDFT
- Electron dynamics capabilities
- Continued development beyond Qbox

### Molecular Dynamics:
- Both BO-MD and CP-MD
- Various thermostats
- Production-quality simulations

## Inputs & Outputs
- **Input formats**:
  - Qball input files (.i)
  - Coordinate files (.sys)
  - Pseudopotential files (.xml)
  - Example inputs in examples/ directory
  
- **Output data types**:
  - Energies and forces
  - MD trajectories
  - Electronic structure data
  - Wavefunction outputs

## Interfaces & Ecosystem
- **Build system**: GNU Autotools (autoconf, automake)
- **Dependencies**: BLAS, LAPACK, ScaLAPACK, FFTW
- **Parallelization**: MPI + OpenMP
- **Platforms**: Linux HPC clusters, Blue Gene systems

## Advanced Features

### RT-TDDFT (Qb@ch):
- Real-time propagation
- Strong-field dynamics
- Continued development at UNC

### HPC Optimization:
- Vendor-specific optimizations (IBM ESSL)
- Custom parallel strategies
- Designed for 10,000+ cores

## Performance Characteristics
- **Speed**: Highly optimized for large core counts
- **Scaling**: Excellent to thousands of MPI ranks
- **Accuracy**: Plane-wave systematic convergence
- **Memory**: Distributed memory model

## Computational Cost
- **Ground state**: Plane-wave standard scaling
- **MD**: Efficient for large trajectories
- **RT-TDDFT**: Time-step dependent
- **Parallelization**: Near-linear scaling on HPC

## Limitations & Known Constraints
- **Compilation**: Requires HPC environment and libraries
- **Learning curve**: HPC expertise helpful
- **Pseudopotentials**: Requires compatible formats
- **Documentation**: Technical focus, less beginner-friendly
- **INQ successor**: LLNL now developing GPU-focused INQ code

## Comparison with Other Codes
- **vs Qbox**: Qb@ll fork with LLNL optimizations and Qb@ch RT-TDDFT
- **vs VASP**: Both plane-wave; Qb@ll open-source, HPC optimized
- **vs QE**: Qb@ll smaller codebase, HPC focus; QE broader features
- **vs SALMON**: Both RT-TDDFT capable; Qb@ll also strong in MD
- **Unique strength**: HPC optimization, combined MD + RT-TDDFT development

## Application Areas
- High-pressure physics
- Warm dense matter
- Extreme conditions simulations
- Materials under shock
- Large-scale molecular dynamics
- Ultrafast dynamics (via Qb@ch)

## Best Practices
- Use example input files as templates
- Configure for target HPC architecture
- Consult LLNL documentation for optimization
- Contact maintainers for support

## Community and Support
- Open-source GPL v3
- LLNL development team
- GitHub repository with 8 contributors
- Contact: Erik Draeger, Xavier Andrade (LLNL)
- Academic development at UNC (Qb@ch)

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/LLNL/qball
2. LLNL Computing: https://computing.llnl.gov/projects/qball
3. Original Qbox: http://qboxcode.org/

**Secondary sources**:
1. Qb@ch development at UNC Chapel Hill
2. LLNL publications using Qb@ll

**Confidence**: VERIFIED
- Repository: ACCESSIBLE (GitHub, LLNL)
- License: GPL v3 (LLNL-CODE-635376)
- Active: 8 contributors, ongoing development

**Verification status**: ✅ VERIFIED
