# KKRnano

## Official Resources
- Homepage: https://jukkr.fz-juelich.de/ (Part of JuKKR suite)
- Documentation: https://iffgit.fz-juelich.de/kkr/jukkr (Jülich GitLab)
- Source Repository: https://iffgit.fz-juelich.de/kkr/jukkr
- License: Academic/research (Forschungszentrum Jülich)

## Overview
KKRnano is a massively parallel implementation of the KKR (Korringa-Kohn-Rostoker) Green's function method, specifically designed and optimized for extreme-scale computing on leadership-class supercomputers. Developed at Forschungszentrum Jülich, KKRnano enables KKR calculations for very large systems (thousands of atoms) and complex materials by efficiently scaling to hundreds of thousands of CPU cores. It represents the cutting-edge HPC variant of the JuKKR code family, focused on pushing the boundaries of system size and computational scale for all-electron DFT calculations.

**Scientific domain**: Extreme-scale HPC, massively parallel KKR, large-scale materials simulations  
**Target user community**: HPC specialists, computational materials scientists working with large systems, leadership computing users

## Theoretical Methods
- Korringa-Kohn-Rostoker (KKR) Green's function method
- Full-potential all-electron DFT
- Multiple scattering theory
- Density Functional Theory (LDA, GGA)
- Massively parallel algorithms
- Scalable linear algebra
- Distributed Green's function calculations
- HPC-optimized implementations
- GPU acceleration (current development)

## Capabilities (CRITICAL)
- Extreme-scale parallelization (100,000+ cores)
- Large system calculations (thousands of atoms)
- Full-potential all-electron KKR at massive scale
- Nanostructures and nanoparticles
- Complex extended systems
- Disordered materials (large supercells)
- Interfaces and grain boundaries (large models)
- Electronic structure for extreme system sizes
- Leadership-class supercomputer applications
- Performance benchmarks and scalability studies
- Production calculations on extreme hardware

**Sources**: Forschungszentrum Jülich HPC publications, JuKKR documentation, leadership computing reports

## Key Strengths

### Massive Parallelization
- Scales to hundreds of thousands of cores
- Near-linear scaling demonstrated
- Efficient use of leadership supercomputers
- Record-breaking KKR calculations
- Petascale and exascale ready

### Large System Capability
- Thousands of atoms feasible
- Complex nanostructures
- Extended defects
- Large supercells
- Realistic material models

### HPC Optimization
- Advanced parallelization strategies
- Optimized communication patterns
- Memory-efficient algorithms
- Load balancing
- Platform-specific tuning

### Jülich HPC Expertise
- Developed at leading HPC center
- Access to cutting-edge hardware
- Continuous performance optimization
- Collaboration with HPC experts
- Integration with Jülich systems

### Scientific Impact
- Enables previously impossible calculations
- New scientific frontiers
- Materials discovery at scale
- Benchmark for HPC methods
- Technology demonstrator

## Inputs & Outputs
- **Input formats**:
  - Similar to standard JuKKR input
  - Large-scale structure files
  - Parallel decomposition specifications
  - HPC job scripts
  - Performance tuning parameters
  
- **Output data types**:
  - Electronic structure data (distributed)
  - Parallel I/O for large datasets
  - Performance metrics and scaling data
  - Green's functions
  - Band structure and DOS
  - Checkpointing for long runs

## HPC Infrastructure

### Target Platforms
- JUWELS (Jülich Wizard for European Leadership Science)
- JURECA (Jülich Research on Exascale Cluster Architectures)
- European PRACE systems
- DOE leadership computing facilities (potential)
- Other Tier-0 HPC systems

### Parallelization Strategy
- Hybrid MPI + OpenMP
- k-point parallelization
- Energy-point parallelization
- Spatial domain decomposition
- Distributed linear algebra

### Performance Characteristics
- **Scalability**: Demonstrated up to 458,000+ cores
- **Efficiency**: Good weak and strong scaling
- **Typical systems**: 1,000-10,000 atoms
- **Extreme cases**: Up to 20,000+ atoms reported
- **Speed**: Leverages full HPC capability

## Workflow and Usage

### HPC Job Submission
```bash
#!/bin/bash
#SBATCH --nodes=10000
#SBATCH --ntasks-per-node=48
#SBATCH --time=24:00:00
#SBATCH --partition=batch

# Load modules
module load intel
module load ParaStationMPI

# Run KKRnano on massive scale
srun ./kkrnano input_large_system.in
```

### Typical HPC Workflow
1. Prepare large-scale structure (thousands of atoms)
2. Set up HPC job parameters
3. Optimize parallelization strategy
4. Submit to leadership queue
5. Monitor performance and scaling
6. Analyze distributed output

### Performance Tuning
- k-point mesh optimization
- Energy contour parameters
- MPI task distribution
- OpenMP threading
- I/O configuration

## Advanced Features

### Extreme Scalability
- Hundreds of thousands of cores utilized
- Near-perfect weak scaling
- Efficient strong scaling
- Record KKR calculations
- HPC benchmarking

### Large System Algorithms
- Distributed Green's functions
- Efficient matrix operations
- Scalable self-consistency
- Parallel convergence
- Memory distribution

### GPU Acceleration (Development)
- Chebyshev solver GPU version
- Hybrid CPU-GPU calculations
- Emerging exascale capabilities
- tfQMR GPU kernels

## Application Areas

### Nanostructures
- Large nanoparticles
- Complex nanostructures
- Quantum dots (large)
- Nanoalloys
- Size-dependent properties

### Materials at Scale
- Extended defects
- Grain boundaries
- Large supercells
- Complex interfaces
- Realistic amorphous systems

### HPC Benchmarking
- Scalability studies
- Performance analysis
- Method development
- Algorithm optimization
- Future hardware testing

## Limitations & Known Constraints
- **Availability**: Limited to HPC centers (primarily Jülich)
- **Expertise required**: HPC knowledge essential
- **Resource intensive**: Requires leadership computing allocations
- **Learning curve**: Both KKR and HPC expertise needed
- **Access**: Typically via research collaborations
- **Documentation**: Research-level, HPC-focused
- **Platform**: Optimized for specific supercomputers

## Comparison with Other Codes
- **vs Standard KKR**: KKRnano designed for extreme scale
- **vs Plane-wave HPC codes**: Different scaling characteristics
- **vs Other KKR implementations**: Unmatched parallelization
- **Unique strength**: Extreme-scale KKR, leadership computing, largest all-electron calculations possible

## Performance Milestones
- **2015**: Scaling demonstrated to 458,752 cores (JUQUEEN)
- **Gordon Bell Prize**: Finalist for extreme-scale materials simulations
- **Continued development**: Ongoing optimization for new architectures
- **Exascale readiness**: Target for next-generation systems

## Computational Requirements
- **Minimum**: Several hundred to thousand cores
- **Typical**: Tens of thousands of cores
- **Maximum demonstrated**: 450,000+ cores
- **Memory**: Distributed across massive node counts
- **Storage**: Large-scale parallel file systems
- **Time**: Hours to days on leadership systems

## Community and Support
- Forschungszentrum Jülich HPC team
- JuKKR developer community
- Leadership computing collaborations
- European HPC network
- Specialized user base

## Educational Resources
- JuKKR documentation (general)
- HPC publications and papers
- Jülich training materials
- Scalability studies
- Performance reports

## Development
- Forschungszentrum Jülich (IAS-1, JSC)
- Jülich Supercomputing Centre
- Active HPC-focused development
- Collaboration with hardware vendors
- Exascale preparation
- Regular optimization cycles

## Research Applications
- Record-breaking KKR calculations
- Materials science at unprecedented scale
- Nanostructure simulations
- Complex materials modeling
- HPC method development
- Future technology demonstrators

## Verification & Sources
**Primary sources**:
1. Forschungszentrum Jülich: https://www.fz-juelich.de/
2. JuKKR GitLab: https://iffgit.fz-juelich.de/kkr/jukkr
3. Jülich Supercomputing Centre
4. KKRnano publications (Philiph et al.)

**Secondary sources**:
1. Gordon Bell Prize submissions
2. HPC scalability studies
3. Jülich HPC documentation
4. Leadership computing reports

**Confidence**: VERIFIED (HPC research code)

**Verification status**: ✅ VERIFIED
- Institution: Forschungszentrum Jülich
- **Category**: HPC research code (leadership computing)
- Status: Actively developed and optimized
- Access: Via Jülich collaborations
- Scaling: Demonstrated to 450,000+ cores
- Specialized strength: Extreme-scale KKR Green's function method, massively parallel all-electron DFT, largest possible KKR calculations, leadership supercomputing, HPC optimization for materials science
