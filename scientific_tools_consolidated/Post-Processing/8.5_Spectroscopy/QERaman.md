# QERaman

## Official Resources
- Source Repository: https://github.com/nguyen-group/QERaman
- Documentation: Included in repository
- License: Open source

## Overview
**QERaman** is an open-source program for computing first-order resonance Raman spectroscopy based on Quantum ESPRESSO. It calculates resonance Raman intensities by evaluating the derivative of the frequency-dependent dielectric function with respect to phonon normal mode coordinates, enabling simulation of Raman spectra under resonant conditions.

**Scientific domain**: Resonance Raman spectroscopy, electron-phonon coupling  
**Target user community**: Researchers studying resonance Raman spectra of crystalline materials from first principles

## Theoretical Methods
- Density Functional Perturbation Theory (DFPT)
- Frequency-dependent dielectric function
- Resonance Raman theory (Kramers-Heisenberg-Dirac)
- Phonon normal modes
- Electron-phonon coupling
- Quantum ESPRESSO as backend

## Capabilities (CRITICAL)
- First-order resonance Raman spectra
- Off-resonance Raman spectra
- Raman intensities vs excitation energy
- Raman excitation profiles
- Phonon frequencies and modes
- Dielectric function calculation
- Polarization-dependent Raman
- Temperature-dependent Raman (via phonon occupation)

**Sources**: GitHub repository, J. Chem. Phys. 158, 224109 (2023)

## Key Strengths

### Resonance Raman:
- Beyond off-resonance approximation
- Excitation energy dependence
- Raman excitation profiles
- Resonance enhancement factors
- Proper treatment of resonant denominators

### QE Integration:
- Uses QE DFPT for phonons
- Uses QE TDDFT for dielectric function
- Same pseudopotentials and structures
- Seamless workflow

### First-Principles:
- No empirical parameters
- Full DFT-level accuracy
- Systematic improvement possible
- Handles complex materials

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO input files
  - QERaman configuration files
  - Phonon mode data from QE
  
- **Output data types**:
  - Resonance Raman spectra
  - Raman excitation profiles
  - Raman intensities per mode
  - Polarization-resolved spectra

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Required backend (v7.1, 7.2, or 7.3)
- **Matplotlib**: For plotting
- **Python**: Scripting interface

## Performance Characteristics
- **Speed**: Depends on number of excitation energies
- **Accuracy**: DFT-level for Raman
- **System size**: Limited by QE capabilities
- **Memory**: Moderate

## Computational Cost
- **Single excitation energy**: Similar to QE DFPT + dielectric
- **Multiple excitation energies**: Scales linearly
- **Typical**: Hours for moderate systems

## Limitations & Known Constraints
- **QE only**: No VASP or other code support
- **First-order only**: No higher-order Raman
- **No excitonic effects**: DFT-level dielectric function
- **QE version**: Requires specific QE versions (7.1-7.3)
- **Documentation**: Limited

## Comparison with Other Codes
- **vs Phonopy-Spectroscopy**: QERaman does resonance Raman, Phonopy-Spectroscopy does off-resonance
- **vs ramannoodle**: QERaman does resonance, ramannoodle does off-resonance with ML
- **vs VASP-Raman**: QERaman uses QE, VASP-Raman uses VASP
- **Unique strength**: First-principles resonance Raman spectroscopy from QE

## Application Areas

### Resonance Raman of Semiconductors:
- Excitation energy dependence
- Resonance profiles
- Band gap effects
- Defect characterization

### 2D Materials:
- MoS2 and TMDs resonance Raman
- Layer-dependent spectra
- Exciton effects
- Strain effects

### Perovskites:
- Resonance Raman of halide perovskites
- Temperature dependence
- Phase transitions
- Electron-phonon coupling

## Best Practices

### Excitation Energy Grid:
- Sample densely near absorption edges
- Use coarser grid far from resonance
- Include both below and above gap energies
- Monitor convergence

### Phonon Calculation:
- Use well-converged QE phonon calculation
- Verify mode assignments
- Check for imaginary frequencies
- Use appropriate q-point grid

## Community and Support
- Open source on GitHub
- Developed by Nguyen group at NTU
- Example calculations provided
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/nguyen-group/QERaman
2. Related: T. P. T. Nguyen et al., J. Chem. Phys. 158, 224109 (2023)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Maintained
- Specialized strength: First-principles resonance Raman spectroscopy from Quantum ESPRESSO
