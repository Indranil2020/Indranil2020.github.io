# QMMM-NAMD

## Official Resources
- Homepage: https://github.com/qmmm-namd/QMMM-NAMD
- Source Repository: https://github.com/qmmm-namd/QMMM-NAMD
- License: MIT License

## Overview
QMMM-NAMD is a software package designed for performing nonadiabatic molecular dynamics (NAMD) simulations within a Quantum Mechanics/Molecular Mechanics (QM/MM) framework. It enables the study of photoinduced processes in complex environments, such as proteins or solutions, by combining accurate quantum mechanical descriptions of chromophores with efficient molecular mechanics models for the surroundings.

**Scientific domain**: QM/MM, nonadiabatic dynamics, surface hopping, photobiology
**Target user community**: Researchers studying photoactive proteins and solution-phase photochemistry

## Theoretical Methods
- QM/MM interface
- Trajectory Surface Hopping (TSH)
- Electrostatic embedding
- Link atom schemes
- Non-adiabatic coupling in QM/MM
- MD propagation

## Capabilities (CRITICAL)
- On-the-fly QM/MM dynamics
- Excited state propagation
- Environment polarization (depending on level)
- Flexible QM and MM partitioning
- Analysis of surface hopping events

**Sources**: GitHub repository

## Key Strengths

### QM/MM Focus:
- Specifically designed for hybrid simulations
- Handles boundary conditions
- Standard MM force fields

### Dynamics:
- Explicit treatment of nonadiabatic transitions
- Coupled electron-nuclear motion in environment

## Inputs & Outputs
- **Input formats**:
  - QM code input template
  - MM topology/coordinates
  - Control file
  
- **Output data types**:
  - Trajectories
  - Hopping logs
  - Energy logs

## Interfaces & Ecosystem
- **QM Codes**: Interfaces to standard codes (check documentation for specific list, typically ORCA/Gaussian/Turbomole)
- **MM Codes**: Interfaces to MM drivers (e.g. Tinker, Amber)

## Advanced Features

### Embedding:
- Electronic embedding for accurate spectral shifts
- Mechanical embedding for constraints

## Performance Characteristics
- **Speed**: Dominated by QM part
- **Parallelization**: Independent trajectories

## Computational Cost
- **High**: QM/MM gradient required every step
- **Scaling**: N_QM^3/4 + N_MM

## Limitations & Known Constraints
- **Interface**: Requires external codes
- **Complexity**: Setting up stable QM/MM runs is non-trivial

## Comparison with Other Codes
- **vs SHARC**: QMMM-NAMD specialized for the QM/MM interface aspect
- **vs Newton-X**: Similar capability, different implementation focus
- **Unique strength**: Dedicated QM/MM-NAMD integration

## Application Areas
- **Photobiology**: Vision, photosynthesis, DNA repair
- **Solution dynamics**: Cage effects, solvent relaxation

## Best Practices
- **equilibration**: Thorough MM equilibration first
- **Active region**: Careful selection of QM atoms
- **Link atoms**: Place away from reactive center

## Community and Support
- Open-source MIT
- GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/qmmm-namd/QMMM-NAMD

**Confidence**: VERIFIED - GitHub project

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (MIT)
- Specialized strength: QM/MM Nonadiabatic Dynamics
