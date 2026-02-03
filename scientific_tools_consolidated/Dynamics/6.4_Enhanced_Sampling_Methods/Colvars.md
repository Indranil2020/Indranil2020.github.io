# Colvars

## Official Resources
- Homepage: https://colvars.github.io/
- Documentation: https://colvars.github.io/colvars-refman-namd/
- Source Repository: https://github.com/Colvars/colvars
- License: LGPL-3.0

## Overview
Colvars (Collective Variables Module) is a software library for molecular simulation that provides a high-performance implementation of sampling algorithms defined on collective variables. It is integrated into NAMD, LAMMPS, GROMACS, and other MD codes.

**Scientific domain**: Collective variables, enhanced sampling, free energy  
**Target user community**: MD users needing collective variable-based sampling

## Theoretical Methods
- Collective variable definitions
- Adaptive biasing force (ABF)
- Metadynamics
- Umbrella sampling
- Steered MD
- Harmonic restraints

## Capabilities (CRITICAL)
- Extensive CV library
- Multiple biasing methods
- Multi-code support (NAMD, LAMMPS, GROMACS)
- Free energy calculations
- Steered MD
- Custom CVs via scripting

## Key Strengths

### CV Library:
- Many predefined CVs
- Combinations and functions
- Custom CVs
- Well-tested

### Multi-code Support:
- NAMD
- LAMMPS
- GROMACS
- VMD
- Consistent interface

## Inputs & Outputs
- **Input formats**:
  - Colvars configuration file
  - MD engine inputs
  
- **Output data types**:
  - CV trajectories
  - Free energy profiles
  - PMF data

## Interfaces & Ecosystem
- **NAMD**: Native integration
- **LAMMPS**: fix colvars
- **GROMACS**: Plugin
- **VMD**: Analysis

## Advanced Features
- **ABF**: Adaptive biasing force
- **Metadynamics**: Well-tempered variant
- **eABF**: Extended ABF
- **Custom CVs**: Scripted variables
- **Multi-walker**: Parallel sampling

## Performance Characteristics
- Efficient C++ implementation
- Low overhead
- Good parallel scaling
- Optimized CVs

## Computational Cost
- CV evaluation fast
- Biasing overhead low
- Depends on method
- Overall: Efficient

## Best Practices
- Choose appropriate CVs
- Validate CV definitions
- Check convergence
- Use appropriate method

## Limitations & Known Constraints
- Configuration syntax learning curve
- Some CVs code-specific
- Documentation spread across codes

## Application Areas
- Protein dynamics
- Ligand binding
- Conformational sampling
- Free energy calculations
- Steered MD

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/Colvars/colvars
2. G. Fiorin et al., Mol. Phys. 111, 3345 (2013)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, LGPL-3.0)
