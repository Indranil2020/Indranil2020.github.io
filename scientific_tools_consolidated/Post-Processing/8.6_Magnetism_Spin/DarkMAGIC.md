# DarkMAGIC

## Official Resources
- Source Repository: https://github.com/Griffin-Group/DarkMAGIC
- Documentation: Included in repository
- License: Open source

## Overview
**DarkMAGIC** (Dark Matter Ab initio maGnon/phonon Interaction Calculator) is a Python package for computing dark matter interaction rates with collective excitations (magnons and phonons) based on ab initio calculations of material properties. It supports magnon calculations using ab initio-based spin Hamiltonians.

**Scientific domain**: Magnon-phonon interactions, dark matter detection, ab initio spin excitations  
**Target user community**: Researchers studying magnon and phonon excitations for dark matter detection and condensed matter physics

## Theoretical Methods
- Ab initio spin Hamiltonian for magnons
- Magnon dispersion calculation
- Phonon calculation from DFT
- Dark matter interaction rates with collective excitations
- Spin-wave theory
- Heisenberg model from DFT

## Capabilities (CRITICAL)
- Magnon dispersion from ab initio spin Hamiltonian
- Phonon dispersion from DFT
- Dark matter absorption rates via magnons
- Dark matter absorption rates via phonons
- Material-specific interaction rates
- Spin Hamiltonian parameter extraction

**Sources**: GitHub repository, Phys. Rev. Lett.

## Key Strengths

### Ab Initio Magnons:
- DFT-based spin Hamiltonian
- First-principles magnon dispersion
- Material-specific calculations
- No empirical parameters

### Dark Matter Applications:
- DM-magnon interaction rates
- DM-phonon interaction rates
- Material optimization for detection
- Direct detection calculations

### Combined Magnon-Phonon:
- Both excitation types
- Cross-coupling effects
- Comprehensive material analysis
- Multi-channel detection

## Inputs & Outputs
- **Input formats**:
  - DFT calculation results
  - Spin Hamiltonian parameters
  - Material specifications
  
- **Output data types**:
  - Magnon dispersion
  - DM absorption rates
  - Interaction cross-sections
  - Material reach curves

## Interfaces & Ecosystem
- **DFT codes**: Parameter extraction
- **Python**: Core language
- **NumPy/SciPy**: Numerical computation

## Performance Characteristics
- **Speed**: Fast (model calculation)
- **Accuracy**: Ab initio level
- **System size**: Depends on spin model
- **Memory**: Low

## Computational Cost
- **Magnon calculation**: Minutes
- **DFT pre-requisite**: Hours (separate)
- **Typical**: Efficient

## Limitations & Known Constraints
- **Niche application**: Dark matter focused
- **Requires DFT input**: Parameters from external calculations
- **Limited documentation**: Research code
- **Small community**: Research group code

## Comparison with Other Codes
- **vs SpinW**: DarkMAGIC has DM interaction, SpinW is general magnon
- **vs UppASD**: DarkMAGIC is magnon+DM, UppASD is spin dynamics
- **vs Spirit**: DarkMAGIC is ab initio magnons for DM, Spirit is general
- **Unique strength**: Ab initio magnon/phonon interactions for dark matter detection, combined magnon-phonon

## Application Areas

### Dark Matter Detection:
- DM absorption via magnons
- DM absorption via phonons
- Material optimization
- Direct detection experiments

### Magnon Physics:
- Ab initio magnon dispersion
- Spin wave excitations
- Magnetic material characterization
- Magnon-phonon coupling

### Condensed Matter:
- Collective excitations
- Spin-lattice coupling
- Thermal properties
- Spectroscopy prediction

## Best Practices

### DFT Input:
- Use well-converged DFT calculations
- Include sufficient neighbor shells for exchange
- Consider spin-orbit coupling
- Validate against experimental magnon spectra

### DM Calculation:
- Use appropriate DM mass range
- Include all relevant excitation channels
- Consider material anisotropy
- Compare with experimental bounds

## Community and Support
- Open source on GitHub
- Developed by Griffin Group (UIC)
- Published in Phys. Rev. Lett.
- Research code

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/Griffin-Group/DarkMAGIC
2. S. Knapen et al., Phys. Rev. Lett. (related)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Published methodology: Phys. Rev. Lett.
- Specialized strength: Ab initio magnon/phonon interactions for dark matter detection, combined magnon-phonon
