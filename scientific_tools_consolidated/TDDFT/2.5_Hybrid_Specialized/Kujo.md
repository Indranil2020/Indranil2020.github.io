# Kujo

## Official Resources
- Homepage: https://github.com/TovCat/Kujo
- Source Repository: https://github.com/TovCat/Kujo
- License: MIT License

## Overview
Kujo is a Python-based software tool designed for analyzing exciton dynamics in organic single crystals. It focuses on calculating exciton couplings and rates using various approximations, enabling the study of singlet fission, triplet fusion, and charge transport in crystalline environments.

**Scientific domain**: Organic crystals, exciton dynamics, singlet fission, charge transport
**Target user community**: Materials scientists studying organic semiconductors

## Theoretical Methods
- Extended Dipole Model (EDM)
- Point Dipole Model
- Transition density coupling
- Marcus theory rates
- Miller-Abrahams rates
- Exciton hopping rates

## Capabilities (CRITICAL)
- Calculation of electronic couplings (J)
- Crystal structure parsing (CIF)
- Neighbor list generation
- Rate calculation for hopping
- Angular dependence analysis
- Switching between coupling models

**Sources**: GitHub repository

## Key Strengths

### Crystal Structure Handling:
- Direct support for periodic crystals
- Automatic neighbor identification
- Supercell generation

### Coupling Models:
- Flexible choice of approximation
- Cutoff-based model switching
- Atomic transition charges

### Python Integration:
- Easy to script
- Integrates with QC output parsing

## Inputs & Outputs
- **Input formats**:
  - CIF crystal files
  - QC output (for transition densities)
  - Configuration file
  
- **Output data types**:
  - Coupling matrices
  - Rate matrices
  - Neighbor lists

## Interfaces & Ecosystem
- **Input**: Gaussian/ORCA/Q-Chem (via cclib/parsing)
- **Language**: Python

## Advanced Features

### Cutoff Handling:
- Distance-dependent model selection
- Optimizes accuracy vs speed

## Performance Characteristics
- **Speed**: Fast (algebraic models)
- **Bottleneck**: Neighbor search in large supercells

## Computational Cost
- **Low**: Post-processing tool
- **QC**: Requires monomer calculations first

## Limitations & Known Constraints
- **Approximation**: Relies on electrostatic models or transition densities
- **Static**: Typical usage is static crystal structure

## Comparison with Other Codes
- **vs VOTCA**: Kujo is lighter, more focused on crystal couplings
- **Unique strength**: Lightweight crystal exciton analysis

## Application Areas
- **Singlet Fission**: Pentacene, rubrene crystals
- **OLEDs**: Host material transport
- **OFETs**: Charge mobility estimates

## Best Practices
- **Monomer QC**: Reliable transition densities
- **Cutoffs**: Test convergence of J with distance
- **Supercell**: Ensure sufficient size for long-range interactions

## Community and Support
- Open-source MIT
- GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/TovCat/Kujo

**Confidence**: VERIFIED - GitHub project

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (MIT)
- Specialized strength: Organic crystal exciton couplings
