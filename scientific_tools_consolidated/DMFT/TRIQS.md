# TRIQS

## Official Resources
- Homepage: https://triqs.github.io/
- Documentation: https://triqs.github.io/triqs/latest/
- Source Repository: https://github.com/TRIQS/triqs
- License: GNU General Public License v3.0

## Overview
TRIQS (Toolbox for Research on Interacting Quantum Systems) is a comprehensive scientific project providing a framework for many-body quantum physics and, in particular, for Dynamical Mean-Field Theory (DMFT) calculations. It consists of C++ libraries with Python interfaces and applications for solving quantum impurity problems and performing DFT+DMFT calculations.

**Scientific domain**: Strongly correlated electron systems, DMFT, many-body physics  
**Target user community**: Condensed matter theorists working on correlated materials

## Theoretical Methods
- Dynamical Mean-Field Theory (DMFT)
- Cluster DMFT extensions
- DFT+DMFT (via DFTTools application)
- Continuous-Time Quantum Monte Carlo (CT-QMC) impurity solvers
- Exact diagonalization solvers
- Hubbard-I approximation
- Green's function formalism
- Self-energy functional theory

## Capabilities (CRITICAL)
- DMFT self-consistency loops
- Quantum impurity solver interfaces (CT-HYB, CT-INT, CT-SEG)
- DFT+DMFT workflows via DFTTools application
- Wannier function downfolding from DFT
- Real-frequency and Matsubara frequency calculations
- Analytical continuation (maximum entropy, Padé)
- Cluster DMFT calculations
- Multi-orbital impurity problems
- Non-local correlations via extended DMFT
- Spectral function calculations
- Python-based workflow scripting
- HDF5-based data storage
- Parallelization support (MPI, OpenMP)

**Sources**: Official TRIQS website, DFTTools documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Python scripts for workflow definition
  - HDF5 archives for Green's functions and self-energies
  - DFT outputs via DFTTools (Wien2k, VASP, Quantum ESPRESSO, ABINIT, Wannier90)
  - Configuration files for impurity solvers
  
- **Output data types**:
  - Green's functions (imaginary time, Matsubara, real frequency)
  - Self-energies
  - Spectral functions
  - Occupations and observables
  - HDF5 archives with full calculation state
  - Python-readable data structures

## Interfaces & Ecosystem
- **TRIQS Applications** (official):
  - TRIQS/cthyb - CT-HYB impurity solver
  - TRIQS/DFTTools - DFT+DMFT interface
  - TRIQS/maxent - Maximum entropy analytical continuation
  - TRIQS/hubbardI - Hubbard-I solver
  - solid_dmft - High-level DFT+DMFT workflows
  
- **DFT code interfaces** (via DFTTools):
  - Wien2k - extensive support
  - VASP - supported via Wannier90 interface
  - Quantum ESPRESSO - supported
  - ABINIT - supported
  - Elk - supported
  - Wannier90 - direct interface for downfolding
  
- **Framework integrations**:
  - AiiDA - workflow automation possible via aiida-triqs (community)
  - Jupyter notebooks - native Python API support
  
- **External solvers**:
  - w2dynamics - can be interfaced
  - Pomerol - exact diagonalization
  - ALPS/CT-HYB - compatibility layer
  
- **Analysis tools**:
  - Python-based post-processing
  - Integration with matplotlib, numpy, scipy
  - HDFView for data inspection

## Limitations & Known Constraints
- **Learning curve**: Steep learning curve; requires understanding of Python, C++, and DMFT theory
- **Installation**: Complex build system with many dependencies; can be challenging to compile
- **Computational cost**: DMFT calculations are inherently expensive; CT-QMC scales poorly with inverse temperature
- **DFT interface setup**: DFT+DMFT requires careful setup of projectors and Wannier functions
- **Memory**: Large memory requirements for multi-orbital problems
- **Real-frequency calculations**: Analytical continuation is ill-posed; results can be unreliable without careful validation
- **Documentation**: While extensive, documentation scattered across main TRIQS and application docs
- **Platform support**: Best supported on Linux; limited Windows support
- **HDF5 version sensitivity**: Can have compatibility issues between different HDF5 versions

## Verification & Sources
**Primary sources**:
1. Official website: https://triqs.github.io/
2. TRIQS documentation: https://triqs.github.io/triqs/latest/
3. DFTTools documentation: https://triqs.github.io/dft_tools/latest/
4. O. Parcollet et al., Comput. Phys. Commun. 196, 398-415 (2015) - TRIQS 1.4
5. P. Seth et al., Comput. Phys. Commun. 200, 274-284 (2016) - TRIQS/cthyb

**Secondary sources**:
1. GitHub repositories: https://github.com/TRIQS
2. TRIQS tutorials and workshop materials
3. solid_dmft documentation for DFT+DMFT workflows
4. Community examples and Jupyter notebooks
5. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (GitHub issues, mailing list, Slack)
- Academic citations: >400 (primary TRIQS papers)
- Ecosystem: Multiple maintained applications
- Workshops: Regular TRIQS schools and tutorials
