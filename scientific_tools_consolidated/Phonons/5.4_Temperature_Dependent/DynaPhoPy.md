# DynaPhoPy

## Official Resources
- Homepage: https://abelcarreras.github.io/DynaPhoPy/
- Source Repository: https://github.com/abelcarreras/DynaPhoPy
- Documentation: https://abelcarreras.github.io/DynaPhoPy/
- License: MIT License

## Overview
DynaPhoPy is a computational code for extracting microscopic anharmonic phonon properties from molecular dynamics simulations using the normal-mode-decomposition technique. It calculates quasiparticle phonon frequencies, linewidths, and lifetimes at finite temperatures.

**Scientific domain**: Anharmonic phonons, temperature-dependent properties, MD analysis  
**Target user community**: Researchers studying temperature-dependent phonon properties

## Theoretical Methods
- Normal mode decomposition
- Velocity autocorrelation analysis
- Spectral energy density
- Quasiparticle phonon frequencies
- Phonon linewidths and lifetimes
- Anharmonic renormalization

## Capabilities (CRITICAL)
- Phonon frequency extraction from MD
- Temperature-dependent frequencies
- Phonon linewidths
- Phonon lifetimes
- Anharmonic effects
- LAMMPS/VASP trajectory support
- Phonopy integration

## Key Strengths

### MD-Based Analysis:
- Direct from trajectories
- Full anharmonicity
- Temperature effects
- No perturbation theory

### Phonopy Integration:
- Uses Phonopy force constants
- Consistent workflow
- Familiar interface
- Well-documented

## Inputs & Outputs
- **Input formats**:
  - LAMMPS trajectories
  - VASP XDATCAR
  - Phonopy force constants
  
- **Output data types**:
  - Phonon frequencies
  - Linewidths
  - Lifetimes
  - Spectral functions

## Interfaces & Ecosystem
- **Phonopy**: Force constants
- **LAMMPS**: MD trajectories
- **VASP**: Ab initio MD
- **Python**: Analysis framework

## Limitations & Known Constraints
- Requires long MD trajectories
- Computational cost
- Resolution limits
- Classical MD effects

## Application Areas
- Anharmonic materials
- High-temperature properties
- Phase transitions
- Thermal transport
- Strongly anharmonic systems

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/abelcarreras/DynaPhoPy
2. A. Carreras et al., Comput. Phys. Commun. 221, 221 (2017)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Documentation: Available
- Academic citations: Well-cited
