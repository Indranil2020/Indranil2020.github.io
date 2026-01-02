# w2dynamics

## Official Resources
- Homepage: https://github.com/w2dynamics/w2dynamics
- Documentation: https://github.com/w2dynamics/w2dynamics/wiki
- Source Repository: https://github.com/w2dynamics/w2dynamics
- License: GNU General Public License v3.0

## Overview
w2dynamics is a continuous-time quantum Monte Carlo (CTQMC) impurity solver for multi-orbital systems within dynamical mean-field theory (DMFT). It provides efficient implementation of hybridization expansion and interaction expansion algorithms with MPI parallelization.

**Scientific domain**: Strongly correlated materials, DMFT calculations, many-body physics  
**Target user community**: Researchers studying strongly correlated electron systems

## Theoretical Methods
- Continuous-time quantum Monte Carlo (CTQMC)
- Hybridization expansion (CT-HYB)
- Interaction expansion (CT-INT)
- Dynamical Mean-Field Theory (DMFT) impurity solver
- Multi-orbital Anderson impurity model
- Worm sampling algorithms
- Maximum entropy analytical continuation
- Legendre polynomial expansion

## Capabilities (CRITICAL)
- DMFT impurity solver for multi-orbital systems
- Single-site and cluster DMFT
- Self-energy calculations
- Green's functions in Matsubara and real frequencies
- Spectral functions via analytical continuation
- Density-density correlations
- Spin and orbital susceptibilities
- Magnetic and orbital order parameters
- General multi-orbital interactions
- Spin-orbit coupling
- Crystal field effects
- Interface to DFT codes for DFT+DMFT
- MPI parallelization

**Sources**: Official w2dynamics documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - HDF5-based input files
  - Configuration files (INI format)
  - Hybridization functions from DMFT loop
  
- **Output data types**:
  - HDF5 output with all observables
  - Self-energies (Matsubara frequencies)
  - Green's functions
  - Spectral functions
  - Occupation numbers and double occupancies
  - Correlation functions

## Interfaces & Ecosystem
- **DMFT frameworks**:
  - solid_dmft - interface for DFT+DMFT
  - TRIQS compatibility (limited)
  - Custom DMFT loops possible
  
- **DFT interfaces**:
  - Works with Wannier functions from Wannier90
  - Interface to VASP, WIEN2k via solid_dmft
  
- **Analysis tools**:
  - Python scripts for post-processing
  - MaxEnt for analytical continuation
  - HDF5 format for data exchange

## Limitations & Known Constraints
- **CTQMC cost**: Expensive; sign problem for some systems
- **Analytical continuation**: MaxEnt introduces uncertainties
- **Multi-orbital complexity**: Many parameters to converge
- **Memory intensive**: Large impurity problems demanding
- **Statistical errors**: Monte Carlo method; error bars on results
- **Learning curve**: DMFT and QMC concepts required
- **Documentation**: Good but assumes DMFT familiarity
- **Platform**: Linux/Unix; MPI required for large calculations

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/w2dynamics/w2dynamics
2. Documentation: https://github.com/w2dynamics/w2dynamics/wiki
3. M. Wallerberger et al., Comput. Phys. Commun. 235, 388 (2019) - w2dynamics code
4. N. Parragh et al., Phys. Rev. B 86, 155158 (2012) - CT-HYB implementation

**Secondary sources**:
1. w2dynamics tutorials and examples
2. Published DFT+DMFT applications
3. solid_dmft documentation
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (GitHub, developers)
- Academic citations: >100
- Active development: Regular updates
