# Socorro

## Official Resources
- **Repository**: https://github.com/sandialabs/socorro (Archived)
- **License**: GNU General Public License (GPL)
- **Status**: Legacy / Archived

## Overview
Socorro is a locally-basis-free, plane-wave Density Functional Theory (DFT) code developed at Sandia National Laboratories. It was engineered for extreme scalability on massively parallel supercomputers. While primary development has ceased, it remains a valuable reference for scalable algorithms, particularly for exact-exchange calculations. It includes Time-Dependent DFT (TDDFT) capabilities for excited states.

**Scientific domain**: Materials science, electronic structure, large-scale systems
**Target user community**: HPC researchers, algorithm developers, legacy project support

## Theoretical Methods
- **Kohn-Sham DFT**: Plane-wave basis set implementation
- **TDDFT**: Time-dependent density functional theory for excitations
- **Exact Exchange**: Novel algorithms for efficient Fock exchange construction
- **Pseudopotentials**: Norm-conserving and Projector Augmented Wave (PAW) support

## Capabilities
- **Massive Scalability**: Demonstrated near-ideal scaling up to ~73,000 cores
- **Hybrid Functionals**: Efficient implementation of exact exchange
- **Spin Polarization**: Collinear spin support
- **Modern Hardware**: Support for OpenMP threading and Intel libraries

## Inputs & Outputs
- **Input formats**: Fortran-style input files/namelists
- **Output data types**:
  - Electronic energies and forces
  - Charge densities
  - Wavefunctions
  - TDDFT spectra (absorption)

## Performance Characteristics
- **Scaling**: Exceptional strong scaling for large systems (e.g., gold supercells)
- **Parallelism**: MPI + OpenMP hybrid
- **Efficiency**: Outperformed commercial codes of its era in core scaling for hybrid functionals

## Computational Cost
- **High Core Count**: Designed to utilize large partitions of supercomputers efficiently.
- **Memory**: Distributed memory model allows for large system sizes.

## Limitations & Known Constraints
- **Development**: Code is archived and no longer actively maintained.
- **Support**: No official user support; relies on legacy documentation.
- **Interface**: Less user-friendly than modern Python-driven codes.

## Comparison with Other Codes
- **vs VASP/QE**: Socorro focuses purely on scalability for specific large problems; VASP/QE are general-purpose with larger ecosystems.
- **vs Conquest**: Both focus on large systems, but Socorro uses standard plane-waves with efficient parallelization rather than O(N) basis sets.

## Citations
- **Primary**: "Socorro: a DFT code for the future?" (Sandia Reports/J. Comp. Phys.)
