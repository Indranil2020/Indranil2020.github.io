# FourPhonon

## Official Resources
- Homepage: https://github.com/FourPhonon/FourPhonon
- Source Repository: https://github.com/FourPhonon/FourPhonon
- Documentation: https://github.com/FourPhonon/FourPhonon/wiki
- License: GPL-3.0

## Overview
FourPhonon is a computational package that extends ShengBTE to calculate four-phonon scattering rates in crystals. It provides exact solutions of the linearized phonon Boltzmann transport equation including four-phonon processes, which are crucial for accurate thermal conductivity predictions in materials with strong anharmonicity.

**Scientific domain**: Thermal transport, four-phonon scattering, lattice thermal conductivity  
**Target user community**: Researchers studying thermal transport in strongly anharmonic materials

## Theoretical Methods
- Four-phonon scattering formalism
- Boltzmann transport equation (BTE)
- Adaptive energy broadening scheme
- Third and fourth-order force constants
- Relaxation time approximation
- Direct solution of linearized BTE

## Capabilities (CRITICAL)
- Four-phonon scattering rate calculations
- Lattice thermal conductivity with 4-phonon processes
- Exact solution of linearized phonon BTE
- Compatible with ShengBTE workflow
- fourthorder.py for 4th-order force constants
- Adaptive broadening for scattering rates
- Temperature-dependent calculations

## Key Strengths

### Four-Phonon Physics:
- Beyond three-phonon approximation
- Essential for strong anharmonicity
- Accurate for high-κ materials
- Captures higher-order effects

### ShengBTE Integration:
- Built on established platform
- Compatible workflow
- Familiar interface
- Proven methodology

### Adaptive Broadening:
- Automatic energy broadening
- Improved numerical stability
- Accurate scattering rates
- Reduced artifacts

## Inputs & Outputs
- **Input formats**:
  - CONTROL file (ShengBTE format)
  - 3rd-order force constants (FORCE_CONSTANTS_3RD)
  - 4th-order force constants (FORCE_CONSTANTS_4TH)
  - Harmonic force constants
  
- **Output data types**:
  - Four-phonon scattering rates
  - Thermal conductivity
  - Mode-resolved properties
  - Relaxation times

## Interfaces & Ecosystem
- **ShengBTE**: Built as extension module
- **thirdorder.py**: 3rd-order force constants
- **fourthorder.py**: 4th-order force constants
- **VASP/QE**: DFT force calculations
- **Phonopy**: Harmonic properties

## Performance Characteristics
- **Computational cost**: Higher than 3-phonon (4th-order scaling)
- **Memory**: Significant for 4th-order tensors
- **Parallelization**: MPI support
- **Accuracy**: Essential for strongly anharmonic systems

## Limitations & Known Constraints
- Higher computational cost than 3-phonon
- Requires 4th-order force constants
- Memory intensive for large systems
- Complex setup for beginners

## Application Areas
- Strongly anharmonic materials
- High thermal conductivity materials (BAs, diamond)
- Thermoelectric materials
- Phase-change materials
- Materials with soft modes

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/FourPhonon/FourPhonon
2. T. Feng et al., Phys. Rev. B 96, 161201(R) (2017)
3. Comput. Phys. Commun. 267, 108033 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-3.0)
- Documentation: Available
- Active development: Yes
- Academic citations: Well-cited
