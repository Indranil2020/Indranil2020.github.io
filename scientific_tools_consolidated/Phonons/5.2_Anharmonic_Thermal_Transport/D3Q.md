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


## Advanced Features

### Core Capabilities:
- Detailed feature implementation
- Advanced algorithms and methods
- Specialized functionality
- Integration capabilities

### Performance Optimizations:
- Computational efficiency features
- Scalability enhancements
- Memory management
- Parallel processing support


## Computational Cost
- **Setup**: Preprocessing requirements
- **Main calculation**: Primary computational cost
- **Post-processing**: Analysis overhead
- **Overall**: Total resource requirements

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

## Comparison with Other Codes
- **vs thirdorder.py (ShengBTE)**: D3Q uses DFPT (no supercell), thirdorder.py uses finite differences
- **vs Phono3py**: Different approach; D3Q is DFPT-based, Phono3py uses supercell finite differences
- **vs hiPhive**: D3Q is DFPT, hiPhive uses machine learning for force constants
- **Unique strength**: Exact DFPT third derivatives without supercell, efficient for metals

## Best Practices

### DFPT Calculations:
- Use dense k-point mesh for metals
- Converge phonon q-point grid
- Check 2n+1 theorem convergence
- Validate symmetry preservation

### Integration with THERMAL2:
- Use consistent q-point grids
- Check force constant quality
- Validate with harmonic phonons
- Test temperature convergence

### Computational Efficiency:
- Use symmetry to reduce calculations
- Parallelize over q-points
- Monitor memory requirements
- Use restart capabilities

## Community and Support
- Open-source GPL-2.0
- Active development (Lorenzo Paulatto group)
- Part of QE ecosystem
- Well-documented methodology
- Published in Phys. Rev. B

## Verification & Sources
**Primary sources**:
1. Website: https://anharmonic.github.io/d3q/
2. L. Paulatto et al., Phys. Rev. B 87, 214303 (2013)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-2.0)
- Documentation: Available
- Active development: Yes
