# QMCPACK

## Official Resources
- Homepage: https://qmcpack.org/
- Documentation: https://qmcpack.readthedocs.io/
- Source Repository: https://github.com/QMCPACK/qmcpack
- License: BSD 3-Clause License

## Overview
QMCPACK is a modern, high-performance quantum Monte Carlo (QMC) code optimized for large-scale electronic structure calculations. It implements variational Monte Carlo (VMC) and diffusion Monte Carlo (DMC) methods for computing ground and excited state properties of molecules and solids with high accuracy. QMCPACK is designed for leadership-class supercomputers and features excellent parallel scaling.

**Scientific domain**: Quantum Monte Carlo, electronic structure, high-accuracy calculations  
**Target user community**: Researchers requiring benchmark-quality electronic structure calculations for molecules and materials

## Theoretical Methods
- Variational Monte Carlo (VMC)
- Diffusion Monte Carlo (DMC)
- Fixed-node DMC
- Reptation quantum Monte Carlo (RQMC)
- Auxiliary-field quantum Monte Carlo (AFQMC)
- Slater-Jastrow trial wavefunctions
- Diffusion Monte Carlo (DMC) with fixed-node approximation
- Reptation Quantum Monte Carlo (RQMC)
- Auxiliary-Field Quantum Monte Carlo (AFQMC)
- Slater-Jastrow wavefunctions
- Multi-determinant expansions
- Backflow transformations
- Twist-averaged boundary conditions
- Finite-size corrections

## Capabilities (CRITICAL)
- Ground-state energies with high accuracy (chemical accuracy achievable)
- Excited-state calculations via DMC with multiple references
- Force calculations for geometry optimization and molecular dynamics
- Periodic boundary conditions for solids
- Open boundary conditions for molecules and clusters
- GPU acceleration (CUDA and HIP)
- Massive parallel scaling (tested up to 100,000+ cores)
- Real-space grids and B-spline orbital representations
- Interface to trial wavefunctions from DFT codes
- Interface to CASINO, PySCF, Quantum ESPRESSO, VASP
- Optimization of Jastrow factors and CI coefficients
- Twist-averaged calculations for finite-size corrections
- Path-integral ground-state calculations
- Response properties via forward walking

**Sources**: Official QMCPACK website, documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - XML input file (qmcpack.in.xml)
  - Trial wavefunctions from DFT (pwscf.h5, ESHDF format)
  - Wavefunctions from quantum chemistry (CASINO, PySCF)
  - Particle set specifications
  - Hamiltonian definitions
  
- **Output data types**:
  - Scalar.dat (energies, variances, acceptance ratios)
  - HDF5 files (wavefunctions, configurations, checkpoints)
  - DMC walker configurations
  - Traces of observables
  - Statistical analysis outputs
  - Optimized wavefunction parameters

## Interfaces & Ecosystem
- **Trial wavefunction sources**:
  - Quantum ESPRESSO - via pw2qmcpack.x converter
  - VASP - supported via converters
  - PySCF - direct Python interface
  - CASINO - wavefunction import
  - Gaussian - via conversion tools
  - GAMESS - supported
  
- **Framework integrations**:
  - Nexus workflow system - official QMCPACK workflow tool
  - ASE - limited integration possible
  - Python interface - for wavefunction preparation
  
- **Post-processing tools**:
  - qmca - analysis of QMCPACK output
  - Nexus analysis utilities
  - Python scripts for data analysis
  
- **HPC optimization**:
  - OpenMP threading
  - MPI parallelization
  - CUDA GPU acceleration (NVIDIA)
  - HIP GPU acceleration (AMD)
  - Hybrid CPU+GPU execution

## Limitations & Known Constraints
- **Fixed-node approximation**: DMC results depend on nodal surface of trial wavefunction; nodal errors cannot be systematically improved
- **Computational cost**: QMC scales as O(N³) to O(N⁴); expensive for large systems
- **Trial wavefunction quality**: Results critically depend on quality of input wavefunction; requires high-quality DFT or CI starting point
- **Pseudopotentials**: Requires specifically designed QMC pseudopotentials; standard DFT pseudopotentials may introduce errors
- **Statistical noise**: QMC results have statistical error bars; long runs needed for high precision
- **Memory**: Large memory requirements for multi-determinant wavefunctions
- **Learning curve**: Steep learning curve; requires understanding of QMC theory and best practices
- **Setup complexity**: Preparing optimal trial wavefunctions and input files is non-trivial
- **System size**: Practical limit ~500-1000 electrons for DMC; smaller for AFQMC

## Verification & Sources
**Primary sources**:
1. Official website: https://qmcpack.org/
2. Documentation: https://qmcpack.readthedocs.io/
3. J. Kim et al., J. Phys.: Condens. Matter 30, 195901 (2018) - QMCPACK code paper
4. P. Kent et al., J. Chem. Phys. 152, 174105 (2020) - QMCPACK AFQMC

**Secondary sources**:
1. QMCPACK tutorials and workshops
2. Nexus documentation for workflow management
3. PySCF integration documentation
4. GPU acceleration papers and benchmarks
5. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (GitHub issues, mailing list)
- Academic citations: >300 (main code papers)
- DOE support: Actively developed and maintained by DOE labs
- GPU support: Verified CUDA and HIP implementations
- Scaling: Demonstrated on leadership-class supercomputers
