# TBPLaS (Tight-Binding Package for Large-scale Simulation)

## Official Resources
- Homepage: https://www.tbplas.net/
- Documentation: https://www.tbplas.net/
- Source Repository: https://github.com/deepmodeling/tbplas
- License: BSD 3-Clause License

## Overview
TBPLaS is a Python package designed for large-scale tight-binding simulations with a focus on performance and scalability. Developed by the DeepModeling community (associated with DeepMD-kit), TBPLaS provides efficient algorithms for electronic structure calculations in tight-binding frameworks, enabling simulations of systems with millions of orbitals. It implements exact diagonalization, tight-binding propagation method (TBPM), and Green's function methods.

**Scientific domain**: Tight-binding, large-scale simulation, HPC, Quantum Transport
**Target user community**: Large-scale TB calculations, HPC users, production simulations

## Theoretical Methods
- Tight-binding Hamiltonian
- Exact Diagonalization (ED)
- Tight-Binding Propagation Method (TBPM)
- Kernel Polynomial Method (KPM)
- Green's function methods
- Sparse matrix techniques
- Band structure & DOS
- Optical conductivity
- Quasi-eigenstates

## Capabilities (CRITICAL)
**Category**: Large-scale Python TB package
- **Scale**: Capable of simulating millions of atoms
- **Methods**: ED, TBPM, KPM, Green's Functions
- **Properties**:
  - Band structure
  - Density of States (DOS)
  - Local DOS (LDOS)
  - Optical conductivity
  - AC/DC conductivity (Kubo formula)
  - Hall conductivity
  - Polarization
- **Performance**:
  - Linear scaling with system size (TBPM)
  - Cython/Fortran extensions
  - MPI/OpenMP parallelization
- **Interface**:
  - Object-oriented Python API
  - Integration with ASE (Atomic Simulation Environment)

**Sources**: Official website, GitHub, arXiv:2209.00806

## Key Strengths

### Large-Scale Simulation:
- Handles systems up to tens of millions of orbitals
- Linear scaling algorithms (TBPM)
- Ideal for disordered systems, quasicrystals, and Moiré superlattices

### Comprehensive Physics:
- Goes beyond simple bands
- Transport properties (conductivity, Hall effect)
- Optical properties
- Spectral functions

### Performance Optimization:
- Critical parts in Cython/Fortran
- Efficient sparse matrix operations
- Optimized for HPC environments

## Status
- **Type**: Large-scale TB package
- **Development**: Active (DeepModeling)
- **Community**: HPC TB users
- **Platform**: Linux, macOS, Windows
- **Integration**: Part of DeepModeling ecosystem (DeepMD-kit, etc.)

## Application Areas
- 2D Materials (Graphene, TMDs)
- Twisted Bilayer Systems (Moiré Physics)
- Disordered Systems (Anderson Localization)
- Quasicrystals
- Quantum Transport in Nanodevices

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.tbplas.net/
2. GitHub: https://github.com/deepmodeling/tbplas
3. Publication: Li, Y., et al. "TBPLaS: a Tight-Binding Package for Large-scale Simulation." arXiv:2209.00806 (2022).

**Confidence**: VERIFIED - Active Research Tool

**Verification status**: ✅ CONFIRMED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- **Note**: Repository moved to `deepmodeling` organization. Confirmed as a high-performance tool for large-scale TB simulations.
