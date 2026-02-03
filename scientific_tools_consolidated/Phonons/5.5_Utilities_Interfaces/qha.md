# qha (Python QHA Package)

## Official Resources
- Homepage: https://github.com/MineralsCloud/qha
- Source Repository: https://github.com/MineralsCloud/qha
- Documentation: https://qha.readthedocs.io/
- Reference: Comput. Phys. Commun. 237, 199 (2019)
- License: GPL-3.0

## Overview
qha is a Python package for quasi-harmonic approximation (QHA) calculations of thermodynamic properties. It computes free energies, thermal expansion, bulk modulus, and other properties as functions of temperature and pressure.

**Scientific domain**: Quasi-harmonic approximation, thermodynamic properties  
**Target user community**: Researchers computing temperature/pressure-dependent properties

## Theoretical Methods
- Quasi-harmonic approximation
- Equation of state fitting
- Free energy minimization
- Thermal expansion
- Grüneisen parameters
- Thermodynamic integration

## Capabilities (CRITICAL)
- Free energy calculations
- Thermal expansion
- Bulk modulus vs T/P
- Heat capacity
- Grüneisen parameters
- Multiple EOS options
- Phonopy integration

## Key Strengths

### QHA Implementation:
- Complete QHA workflow
- Multiple volumes
- Temperature/pressure
- Thermodynamic properties

### Phonopy Compatible:
- Uses Phonopy output
- Familiar workflow
- Well-documented
- Active development

## Inputs & Outputs
- **Input formats**:
  - Phonopy thermal properties
  - Volume-energy data
  - Configuration files
  
- **Output data types**:
  - Free energies
  - Thermal expansion
  - Bulk modulus
  - Heat capacity
  - Grüneisen parameters

## Interfaces & Ecosystem
- **Phonopy**: Phonon input
- **Python**: Analysis framework
- **matplotlib**: Plotting

## Limitations & Known Constraints
- QHA approximation limits
- Requires multiple volumes
- Computational cost
- Anharmonic effects limited

## Application Areas
- Thermodynamic properties
- Phase diagrams
- High-pressure studies
- Thermal expansion
- Geophysics applications

## Comparison with Other Codes
- **vs Phonopy native QHA**: qha provides more EOS options and analysis
- **vs Gibbs2**: qha is Python-based, more flexible
- **vs thermo_pw**: qha works with any phonon code via Phonopy
- **Unique strength**: Comprehensive QHA with multiple EOS models

## Best Practices

### Volume Sampling:
- Use at least 5-7 volume points
- Cover ±5-10% volume range
- Ensure smooth E-V curve
- Check for imaginary modes

### Phonon Calculations:
- Use consistent settings across volumes
- Converge q-point mesh
- Check for negative frequencies
- Validate harmonic approximation

### EOS Fitting:
- Try multiple EOS models
- Check fitting quality
- Validate pressure range
- Compare with experiments

### Temperature Range:
- Start from low temperature
- Extend to relevant T range
- Check QHA validity limits
- Monitor anharmonic effects

## Community and Support
- Open-source GPL-3.0
- Active development (MineralsCloud)
- Published methodology (CPC 2019)
- ReadTheDocs documentation
- Examples and tutorials

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/MineralsCloud/qha
2. T. Qin et al., Comput. Phys. Commun. 237, 199 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-3.0)
- Published paper
- Active development
