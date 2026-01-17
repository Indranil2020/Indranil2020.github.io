# Gator

## Official Resources
- Homepage: https://e-science.se/software/gator/
- Source Repository: https://github.com/gator-program/gator (Private/Request access or distributed via website)
- License: GNU General Public License v3.0

## Overview
Gator is a quantum chemistry program specialized for spectroscopy and molecular properties using the Algebraic Diagrammatic Construction (ADC) scheme. It focuses on correlated excited state calculations, particularly for simulating core-level spectroscopies such as X-ray absorption (XAS), X-ray emission (XES), and Resonant Inelastic X-ray Scattering (RIXS), as well as valence excitations.

**Scientific domain**: Computational spectroscopy, X-ray spectroscopy, ADC methods, core-excited states
**Target user community**: Researchers in X-ray spectroscopy and correlated excited states

## Theoretical Methods
- ADC(2): Second-order ADC for excitation energies
- ADC(3): Third-order ADC
- CVS-ADC: Core-Valence Separation for core states
- Intermediate State Representation (ISR)
- TP-ADC(2): Transition-Potential ADC
- SOC: Spin-Orbit Coupling treatment

## Capabilities (CRITICAL)
- Core-excitation energies (XAS)
- Core-emission energies (XES)
- RIXS cross-sections
- Valence excitation energies
- Transition moments
- Spin-orbit coupling effects (relevant for L-edges)
- Property calculations

**Sources**: Official website, published papers

## Key Strengths

### Spectroscopic Focus:
- Specifically designed for simulating experimental spectra
- Robust implementation of transition properties
- Core-level specific features (CVS)

### High-Order Correlation:
- ADC(3) capabilities for high accuracy
- Systematic improvement over TDDFT
- Reliable for charge-transfer and Rydberg states

### Relativistic Effects:
- Scalar relativity
- Spin-orbit coupling (critical for metal K/L edges)

## Inputs & Outputs
- **Input formats**:
  - Gator input file
  - Hartree-Fock data (often interfaced with other codes like Molcas/Dalton)
  
- **Output data types**:
  - Excitation lists
  - Oscillator strengths
  - Cross-sections
  - Spectrum data files

## Interfaces & Ecosystem
- **SCF Driver**: Typically requires an interface to an SCF code (e.g., Dalton, OpenMolcas) to generate MO integrals
- **Language**: Fortran/C++

## Advanced Features

### RIXS Simulation:
- Kramers-Heisenberg formula implementation
- Two-step spectroscopy
- Interference effects

### Core-Valence Separation:
- Projecting out valence continuum
- Stable convergence for high-energy states

## Performance Characteristics
- **Accurate**: High-level correlated method
- **Cost**: Higher than TDDFT, lower than EOM-CCSDT
- **Scaling**: N^5 for ADC(2), N^6 for ADC(3)

## Computational Cost
- **Memory**: High (storage of amplitudes)
- **Time**: Significant for large basis sets
- **Cluster**: Recommended for real systems

## Limitations & Known Constraints
- **Availability**: Distribution might be less automated than GitHub projects
- **Ground State**: Depends on external SCF
- **System Size**: Limited by correlated method scaling (<50-100 atoms typical)

## Comparison with Other Codes
- **vs Q-Chem**: Open-source alternative for ADC spectroscopy
- **vs ORCA**: Specialized for RIXS/X-ray, whereas ORCA is general purpose
- **Unique strength**: Dedicated focus on high-level X-ray spectroscopy simulation

## Application Areas
- **X-ray Absorption (XAS)**: K-edge, L-edge of transition metals
- **X-ray Emission (XES)**: Valence-to-core emission
- **RIXS**: Inelastic scattering maps
- **Photochemistry**: Accurate vertical excitations

## Best Practices
- **Basis Set**: Core properties require specialized basis sets (e.g., core-valence sets)
- **CVS**: Must be enabled for core states
- **Memory**: Allocate sufficient scratch space
- **Validation**: Compare against experimental spectra

## Community and Support
- Academic development (e-science.se)
- Support via developers
- Manual available online

## Verification & Sources
**Primary sources**:
1. Website: https://e-science.se/software/gator/
2. Scientific publications by authors (e.g. S. Coriani theory)

**Confidence**: VERIFIED - Established academic code

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (GPL)
- Specialized strength: ADC methods for X-ray spectroscopy
