# ThermoPW

## Official Resources
- Homepage: https://dalcorso.github.io/thermo_pw/
- Documentation: https://dalcorso.github.io/thermo_pw/thermo_pw_help.html
- Source Repository: https://github.com/dalcorso/thermo_pw
- License: GNU General Public License v2 (part of Quantum ESPRESSO)

## Overview
**ThermoPW** is a Fortran driver for Quantum ESPRESSO that enables parallel and automated computation of material thermodynamic and spectroscopic properties. It leverages QE routines for calculating dielectric properties, infrared spectra, Raman spectra, elastic constants, and thermodynamic quantities in a streamlined, high-throughput workflow.

**Scientific domain**: Thermodynamic properties, vibrational spectroscopy, dielectric properties, elastic properties  
**Target user community**: Researchers needing automated high-throughput calculation of material properties from DFT

## Theoretical Methods
- Density Functional Perturbation Theory (DFPT)
- Density Functional Theory (DFT)
- Phonon calculations at Γ-point and throughout BZ
- Dielectric tensor calculation
- Born effective charges
- Infrared intensities
- Raman tensors
- Elastic constants
- Quasi-harmonic approximation (QHA)
- Thermal expansion

## Capabilities (CRITICAL)
- Infrared (IR) spectra
- Raman spectra (non-resonant)
- Dielectric function (real and imaginary parts)
- Loss function
- Refractive index
- Born effective charges
- Elastic constants and compliance
- Phonon dispersion at Γ
- Thermodynamic properties (free energy, entropy, heat capacity)
- Thermal expansion (QHA)
- Grüneisen parameters
- High-throughput mode
- Parallel execution over multiple calculations

**Sources**: Official documentation, Comput. Phys. Commun. 293, 108905 (2023)

## Key Strengths

### Automated Workflow:
- Single input for multiple properties
- Automatic convergence testing
- Parallel execution over k-points and strains
- High-throughput capable
- No manual step-by-step execution

### Comprehensive Property Coverage:
- IR and Raman spectra
- Dielectric properties
- Elastic properties
- Thermodynamic properties
- All from same calculation set

### QE Integration:
- Uses QE's well-tested DFPT routines
- Compatible with all QE functionals
- Same pseudopotentials and basis
- Seamless integration

### Parallelization:
- Parallel over images (different strains/volumes)
- Parallel over k-points
- Efficient resource usage
- Scalable to HPC

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO input files
  - Thermo_pw specific namelists
  - Structure files
  
- **Output data types**:
  - IR spectra (frequencies and intensities)
  - Raman spectra (frequencies and cross-sections)
  - Dielectric function vs frequency
  - Elastic constants tensor
  - Thermodynamic quantities vs temperature
  - Phonon frequencies at Γ

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Core engine
- **Phonopy**: Can use thermo_pw as alternative
- **Standard plotting**: Text output for gnuplot/python

## Performance Characteristics
- **Speed**: Efficient parallel execution
- **Accuracy**: Same as QE DFPT
- **System size**: Limited by QE capabilities
- **Parallelization**: Excellent (image + k-point parallel)

## Computational Cost
- **IR/Raman**: Same as QE DFPT phonon calculation
- **Elastic**: 6 strain calculations
- **QHA**: Multiple volume calculations
- **Typical**: Hours for moderate systems

## Limitations & Known Constraints
- **QE only**: No support for other DFT codes
- **Γ-point phonons**: No full BZ dispersion (use Phonopy for that)
- **Non-resonant Raman**: No resonance effects
- **No BSE**: No excitonic effects in optical spectra
- **Learning curve**: Additional namelists beyond QE

## Comparison with Other Codes
- **vs Phonopy-Spectroscopy**: ThermoPW is integrated with QE, Phonopy-Spectroscopy is post-processing
- **vs VASP Raman**: ThermoPW uses QE, VASP-Raman uses VASP
- **vs epsilon/QE**: ThermoPW automates what epsilon.x does manually
- **Unique strength**: Automated high-throughput IR/Raman/dielectric/elastic/thermodynamic properties from QE

## Application Areas

### Vibrational Spectroscopy:
- IR spectra of crystals
- Raman spectra of crystals
- Mode assignments
- Temperature-dependent spectra

### Dielectric Properties:
- Static and frequency-dependent dielectric constant
- Refractive index
- Loss function
- Polarizability

### Elastic Properties:
- Elastic constants
- Bulk and shear moduli
- Anisotropy
- Sound velocities

### Thermodynamics:
- Free energy vs temperature
- Heat capacity
- Thermal expansion
- Phase stability

## Best Practices

### Convergence:
- Test k-point convergence for dielectric
- Use appropriate smearing for metals
- Validate phonon frequencies
- Check Raman tensor convergence

### Parallelization:
- Use image parallelism for elastic/QHA
- Use k-point parallelism for individual calculations
- Balance resources between images and k-points

### Output Analysis:
- Use provided plotting scripts
- Compare IR/Raman with experiment
- Validate dielectric constants
- Cross-check elastic constants

## Community and Support
- Open source (GPL v2, part of QE ecosystem)
- Active development by A. Dal Corso
- Documentation and tutorials available
- Part of Quantum ESPRESSO distribution

## Verification & Sources
**Primary sources**:
1. Official website: https://dalcorso.github.io/thermo_pw/
2. A. Dal Corso, Comput. Phys. Commun. 293, 108905 (2023)
3. GitHub repository: https://github.com/dalcorso/thermo_pw

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub, GPL v2)
- Community support: Active (QE community)
- Academic citations: Growing
- Active development: Ongoing
- Specialized strength: Automated high-throughput IR/Raman/dielectric/elastic/thermodynamic properties from QE
