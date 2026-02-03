# THERMAL2

## Official Resources
- Homepage: https://anharmonic.github.io/thermal2/
- Source Repository: https://github.com/anharmonic/thermal2
- Documentation: https://anharmonic.github.io/thermal2/
- License: GPL-2.0

## Overview
THERMAL2 is a suite of codes for computing lattice thermal conductivity and related anharmonic properties from first principles. It works with third-order force constants from D3Q or other sources to solve the phonon Boltzmann transport equation.

**Scientific domain**: Thermal transport, phonon BTE, lattice thermal conductivity  
**Target user community**: Researchers computing thermal conductivity from first principles

## Theoretical Methods
- Phonon Boltzmann Transport Equation
- Relaxation time approximation
- Variational solution
- Wigner transport equation
- Three-phonon scattering
- Isotope scattering

## Capabilities (CRITICAL)
- Lattice thermal conductivity
- Phonon lifetimes and linewidths
- Variational BTE solution
- Wigner conductivity (quantum corrections)
- Isotope scattering
- Grain boundary scattering
- Temperature-dependent properties

## Key Strengths

### Multiple Methods:
- RTA and variational
- Wigner corrections
- Various scattering mechanisms
- Flexible approach

### D3Q Integration:
- Seamless workflow
- QE compatibility
- DFPT force constants
- Consistent methodology

## Inputs & Outputs
- **Input formats**:
  - D3Q force constants
  - Dynamical matrices
  - Configuration files
  
- **Output data types**:
  - Thermal conductivity tensor
  - Phonon lifetimes
  - Scattering rates
  - Mode contributions

## Interfaces & Ecosystem
- **D3Q**: Third-order force constants
- **Quantum ESPRESSO**: DFT/DFPT
- **q2r.x**: Force constant processing

## Limitations & Known Constraints
- Primarily for D3Q workflow
- Requires anharmonic force constants
- Complex for beginners
- QE-centric

## Application Areas
- Thermal conductivity predictions
- Thermoelectric materials
- Phonon engineering
- Heat management materials

## Verification & Sources
**Primary sources**:
1. Website: https://anharmonic.github.io/thermal2/
2. L. Paulatto et al., Phys. Rev. B 91, 054304 (2015)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-2.0)
- Documentation: Available
- Active development: Yes
