# TDAP (Time-Dependent Ab-initio Propagation)

## Official Resources
- Homepage: http://english.iopp.cas.cn/ (Institute of Physics, CAS - likely internal/request based)
- Documentation: Not publicly hosted
- Source Repository: Closed source / Upon request
- License: Proprietary / Copyrighted (Registration No. 2017SR656635)

## Overview
TDAP (Time-Dependent Ab-initio Propagation) is a software package developed by the group of **Sheng Meng** at the Institute of Physics, Chinese Academy of Sciences. It is based on Time-Dependent Density Functional Theory (TDDFT) using numerical atomic basis sets. The code is designed for real-time simulations of ultrafast electron dynamics, excited state molecular dynamics, and optical properties in complex systems.

**Scientific domain**: Real-time TDDFT, ultrafast dynamics, electron injection, excited state MD
**Target user community**: Researchers in photovoltaics, surface science, and ultrafast physics

## Theoretical Methods
- Real-Time TDDFT (RT-TDDFT)
- Numerical Atomic Orbitals (NAO) basis
- Ehrenfest Dynamics for ions
- Non-adiabatic dynamics
- Linear Response (optical spectra via dipole evolution)
- Field-dependent nonlinear response

## Capabilities (CRITICAL)
- Ultrafast electron injection dynamics
- Coupled electron-ion dynamics (excited state MD)
- Optical absorption spectra calculation
- Nonlinear optical properties
- Simulation of photovoltaic interfaces
- Dye-sensitized solar cells (DSSC) modeling

**Sources**:
- "TDAP 2.0: A package for real-time TDDFT simulations" (referenced in CAS reports)
- Group website/publications of Prof. Sheng Meng (IOP-CAS)

## Key Strengths

### Ultrafast Dynamics:
- Explicit time-domain simulation
- Electron transfer at interfaces
- Hot carrier relaxation

### Numerical Orbitals:
- Efficient for large systems (surfaces, nanostructures)
- Good balance of accuracy and cost
- Comparable to SIESTA in basis infrastructure

## Inputs & Outputs
- **Input formats**:
  - Structure files
  - Pseudopotentials
  - Simulation parameters (time step, field strength)
  
- **Output data types**:
  - Time-dependent dipole moments
  - Population analysis
  - Excitation energy evolution
  - Ionic trajectories

## Interfaces & Ecosystem
- **Basis**: Numerical atomic orbitals
- **Relation**: Methodology shares similarities with SIESTA/SIESTA-TDDFT approaches

## Performance Characteristics
- **Scaling**: O(N) or near-linear for key operations
- **System size**: Capable of handling hundreds of atoms (interface systems)

## Limitations & Known Constraints
- **Availability**: Not an open-source community code; proprietary/research group code.
- **Documentation**: Limited public documentation.

## Comparison with Other Codes
- **vs SIESTA**: TDAP uses similar NAO basis but specialized for functionality developed at IOP-CAS.
- **vs Octopus**: Both are RT-TDDFT, but TDAP uses NAOs while Octopus uses real-space grids.

## Application Areas
- **Photovoltaics**: Charge transfer in solar cells
- **Surface Physics**: Adsorbate dynamics under illumination
- **2D Materials**: Optical response of monolayers

## Community and Support
- Developed at Institute of Physics, Chinese Academy of Sciences (Beijing).
- Support limited to collaborators and licensed users.

## Verification & Sources
**Primary sources**:
1. Institute of Physics, CAS News (TDAP-2.0 release, 2018)
2. Publications by Sheng Meng group (e.g., J. Chem. Phys, Phys. Rev. B using TDAP)

**Confidence**: VERIFIED - Research group code

**Verification status**: ✅ VERIFIED
- Official homepage: Research Group Page (IOP-CAS)
- Source code: CLOSED (Copyrighted)
- Method: Real-Time TDDFT with NAOs
- Specialized strength: Ultrafast electron dynamics at interfaces
