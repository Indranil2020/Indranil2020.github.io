# D3Q

## Official Resources
- Homepage: https://anharmonic.github.io/d3q/
- Source Repository: https://github.com/anharmonic/d3q
- Documentation: https://anharmonic.github.io/d3q/
- License: GPL-2.0

## Overview
D3Q is a code for computing third-order anharmonic force constants from density functional perturbation theory (DFPT) within the Quantum ESPRESSO framework. It enables calculation of phonon-phonon scattering rates and thermal transport properties.

**Scientific domain**: Anharmonic phonons, third-order force constants, DFPT  
**Target user community**: Quantum ESPRESSO users studying anharmonic properties

## Theoretical Methods
- Density Functional Perturbation Theory (DFPT)
- Third-order force constants (2n+1 theorem)
- Phonon-phonon scattering
- Anharmonic perturbation theory
- Wavefunction perturbation recomputation

## Capabilities (CRITICAL)
- Third-order force constant calculation
- DFPT-based anharmonic properties
- Integration with QE ph.x
- Phonon linewidth calculations
- Thermal conductivity (via THERMAL2)
- Efficient 2n+1 implementation

## Key Strengths

### DFPT Approach:
- Exact third derivatives
- No supercell needed
- Efficient for metals
- Systematic accuracy

### QE Integration:
- Works with ph.x
- Familiar workflow
- Well-tested
- Active development

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO files
  - ph.x dynamical matrices
  - pw.x wavefunctions
  
- **Output data types**:
  - Third-order force constants
  - Anharmonic matrices
  - Input for THERMAL2

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Primary integration
- **THERMAL2**: Thermal conductivity
- **ph.x**: Phonon calculations
- **pw.x**: Ground state

## Limitations & Known Constraints
- QE-specific
- Requires DFPT expertise
- Memory intensive
- Complex setup

## Application Areas
- Anharmonic phonon properties
- Thermal conductivity calculations
- Phonon linewidths
- Materials with strong anharmonicity

## Verification & Sources
**Primary sources**:
1. Website: https://anharmonic.github.io/d3q/
2. L. Paulatto et al., Phys. Rev. B 87, 214303 (2013)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-2.0)
- Documentation: Available
- Active development: Yes
