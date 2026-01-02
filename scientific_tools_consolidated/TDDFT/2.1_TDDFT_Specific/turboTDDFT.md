# turboTDDFT (Quantum ESPRESSO TDDFT)

## Official Resources
- Homepage: https://www.quantum-espresso.org/
- Documentation: https://github.com/dceresoli/qe-gipaw (Legacy turboTDDFT)
- Source Repository: Part of Quantum ESPRESSO
- License: GNU GPL

## Overview
turboTDDFT is a legacy TDDFT module for Quantum ESPRESSO that implements time-dependent density functional theory for calculating optical absorption spectra, excitation energies, and dynamic polarizabilities. Historically part of the QE ecosystem, turboTDDFT functionality has largely been superseded by newer implementations (turboEELS, turbo_spectrum.x) within Quantum ESPRESSO. It uses linear-response TDDFT with plane-wave basis sets for periodic systems.

**Scientific domain**: TDDFT, optical properties, linear response  
**Target user community**: Quantum ESPRESSO users, solid-state spectroscopists

## Theoretical Methods
- Time-Dependent Density Functional Theory (TDDFT)
- Linear-response formalism
- Plane-wave basis sets
- Pseudopotentials
- LDA and GGA functionals
- Periodic systems
- Casida equation
- Optical absorption

## Capabilities (CRITICAL)
**Note**: Legacy module, functionality now in newer QE tools.
- Linear-response TDDFT
- Optical absorption spectra
- Excitation energies
- Oscillator strengths
- Dynamic polarizability
- Plane-wave implementation
- Periodic systems
- Integration with Quantum ESPRESSO

**Sources**: Quantum ESPRESSO documentation

## Key Characteristics

### Quantum ESPRESSO Integration:
- Part of QE ecosystem
- Uses QE wavefunctions
- Plane-wave basis
- Pseudopotential framework
- Standard QE workflow

### Linear-Response TDDFT:
- Casida-like approach
- Excitation energies
- Absorption spectra
- Standard formalism
- Production calculations

### Legacy Status:
- Historical module
- Superseded by turbo_spectrum.x
- turboEELS for EELS
- Newer implementations preferred
- Documentation limited

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO input
  - Ground-state wavefunctions
  - TDDFT parameters
  
- **Output data types**:
  - Excitation energies
  - Oscillator strengths
  - Absorption spectra
  - Polarizabilities

## Interfaces & Ecosystem
- **Quantum ESPRESSO**:
  - Integrated module
  - Uses pw.x output
  - QE input format
  - Standard workflow
  
- **Superseded By**:
  - turbo_spectrum.x
  - turboEELS
  - turbo_lanczos.x
  - Newer QE TDDFT tools

## Workflow

### Historical Workflow:
1. Run Quantum ESPRESSO SCF
2. Generate ground-state wavefunctions
3. Run turboTDDFT
4. Analyze optical spectra

### Modern Alternative:
- Use turbo_spectrum.x or turboEELS in current QE versions

## Legacy Code Status

### Historical Context:
- Early QE TDDFT implementation
- Pioneering plane-wave TDDFT
- Important historical tool
- Now superseded

### Current Recommendations:
- Use turbo_spectrum.x for optical spectra
- Use turboEELS for electron energy loss
- Use turbo_lanczos.x for large systems
- Consult QE documentation for current tools

## Limitations & Known Constraints
- **Legacy status**: Superseded by newer tools
- **Documentation**: Limited for old version
- **Maintenance**: No longer actively developed
- **Support**: Community support limited
- **Recommendation**: Use newer QE TDDFT modules

## Modern QE TDDFT Tools

### turbo_spectrum.x:
- Current optical absorption tool
- Improved algorithms
- Better performance
- Active development

### turboEELS:
- Electron energy loss spectroscopy
- Modern implementation
- Comprehensive features
- Well-documented

### turbo_lanczos.x:
- Large system TDDFT
- Efficient algorithms
- Scalable
- Production quality

## Quantum ESPRESSO TDDFT Ecosystem

### Current Tools:
- **turbo_spectrum.x**: Optical absorption
- **turboEELS**: EELS calculations
- **turbo_lanczos.x**: Large-scale TDDFT
- **turbo_davidson.x**: Alternative solver

### Recommendation:
For Quantum ESPRESSO TDDFT calculations, users should consult the current QE documentation and use the actively maintained TDDFT modules rather than legacy turboTDDFT.

## Historical Significance
- Early plane-wave TDDFT
- QE ecosystem pioneer
- Enabled TDDFT in QE
- Foundation for current tools
- Important historical contribution

## Verification & Sources
**Primary sources**:
1. Quantum ESPRESSO: https://www.quantum-espresso.org/
2. QE documentation (historical references)
3. Legacy turboTDDFT code repositories

**Secondary sources**:
1. TDDFT literature
2. Quantum ESPRESSO user community
3. Historical QE publications

**Confidence**: VERIFIED - Legacy module

**Verification status**: ✅ VERIFIED (Legacy)
- Status: **LEGACY MODULE** - Superseded by newer QE TDDFT tools
- Quantum ESPRESSO: CONFIRMED
- Current recommendation: Use turbo_spectrum.x, turboEELS, or turbo_lanczos.x
- Documentation: Consult current Quantum ESPRESSO documentation
- **For TDDFT in Quantum ESPRESSO**: Use actively maintained modules in current QE versions
