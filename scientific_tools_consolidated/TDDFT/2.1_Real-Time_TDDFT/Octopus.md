# Octopus

## Official Resources
- Homepage: https://octopus-code.org/
- Documentation: https://octopus-code.org/documentation/
- Source Repository: https://gitlab.com/octopus-code/octopus
- License: GNU General Public License v2.0

## Overview
Octopus is a scientific program for the ab initio simulation of electron-ion dynamics using time-dependent density-functional theory (TDDFT) and real-space grids. It specializes in real-time propagation for studying ultrafast processes, optical properties, and electron dynamics.

**Scientific domain**: Ultrafast dynamics, optical properties, electron-ion dynamics, strong-field physics  
**Target user community**: Researchers studying time-dependent phenomena, optical response, laser-matter interaction

## Theoretical Methods
- Time-Dependent Density Functional Theory (TDDFT)
- Real-time TDDFT propagation
- Density Functional Theory (DFT) for ground state
- Real-space finite-differences method
- LDA, GGA, hybrid functionals
- Optimal control theory (OCT)
- Ehrenfest molecular dynamics
- Multi-system calculations (electron-ion, photon-electron)
- Non-adiabatic dynamics
- Casida equation for linear response

## Capabilities (CRITICAL)
- Ground-state DFT calculations
- Real-time TDDFT for electron dynamics
- Optical absorption spectra
- Time-resolved spectroscopy simulation
- Strong-field physics (high-harmonic generation)
- Photoionization and photoemission
- Ehrenfest molecular dynamics
- Optimal control for laser pulse design
- Non-linear optical properties
- Plasmonic excitations
- Photon-electron coupling
- Multi-component systems
- Geometry optimization
- Vibrational analysis
- GPU acceleration for real-time propagation

**Sources**: Official Octopus documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - inp file (Octopus input format)
  - XYZ coordinate files
  - Pseudopotential files
  
- **Output data types**:
  - Time-dependent observables
  - Absorption spectra
  - Electronic densities (real-space grids)
  - Time-propagated wavefunctions
  - Ehrenfest trajectories
  - Photoelectron spectra

## Interfaces & Ecosystem
- **Framework integrations**:
  - LibXC for exchange-correlation functionals
  - ETSF I/O for file formats
  
- **Visualization**:
  - Compatible with standard grid visualization
  - Custom analysis scripts
  
- **GPU support**:
  - CUDA and OpenCL acceleration
  - Significant speedup for propagation

## Limitations & Known Constraints
- **Real-space grids**: Require convergence testing for grid spacing
- **Pseudopotentials**: Limited to norm-conserving
- **System size**: Real-time TDDFT expensive; ~100-500 atoms typical
- **Time propagation**: Long simulations memory and time intensive
- **k-point sampling**: Best for finite systems or Gamma-point
- **Hybrid functionals**: Computationally expensive
- **Learning curve**: TDDFT concepts require understanding
- **Parallelization**: MPI and GPU but efficiency varies
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://octopus-code.org/
2. Documentation: https://octopus-code.org/documentation/
3. GitLab repository: https://gitlab.com/octopus-code/octopus
4. A. Castro et al., Phys. Status Solidi B 243, 2465 (2006) - Octopus code
5. X. Andrade et al., J. Phys.: Condens. Matter 24, 233202 (2012) - Real-space TDDFT

**Secondary sources**:
1. Octopus tutorials and workshops
2. Published ultrafast dynamics applications
3. Strong-field physics studies
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitLab)
- Community support: Active (mailing list, GitLab)
- Academic citations: >500 (main papers)
- Active development: Regular releases, GPU optimization
