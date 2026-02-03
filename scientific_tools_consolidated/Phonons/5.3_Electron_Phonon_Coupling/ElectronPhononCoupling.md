# ElectronPhononCoupling (EPC)

## Official Resources
- Homepage: https://github.com/GkAntonius/ElectronPhononCoupling
- Source Repository: https://github.com/GkAntonius/ElectronPhononCoupling
- License: GPL-3.0

## Overview
ElectronPhononCoupling is a Python module for analyzing electron-phonon coupling quantities computed with Abinit. It provides tools for computing temperature-dependent band structures, zero-point renormalization, and other electron-phonon related properties.

**Scientific domain**: Electron-phonon coupling, temperature-dependent electronic structure  
**Target user community**: Abinit users studying electron-phonon interactions

## Theoretical Methods
- Allen-Heine-Cardona theory
- Temperature-dependent band gaps
- Zero-point renormalization (ZPR)
- Fan-Migdal self-energy
- Debye-Waller contributions
- Electron-phonon matrix elements

## Capabilities (CRITICAL)
- Temperature-dependent band structures
- Zero-point renormalization
- Band gap temperature dependence
- Electron-phonon self-energy
- Fan and Debye-Waller terms
- Spectral functions
- Integration with Abinit output

## Key Strengths

### Abinit Integration:
- Direct use of Abinit output
- DFPT electron-phonon data
- Consistent methodology
- Well-tested workflow

### Temperature Effects:
- Full temperature dependence
- Zero-point motion
- Quantum effects
- Accurate predictions

## Inputs & Outputs
- **Input formats**:
  - Abinit netCDF files
  - Electron-phonon matrix elements
  - Phonon frequencies
  
- **Output data types**:
  - Temperature-dependent bands
  - Renormalized gaps
  - Self-energies
  - Spectral functions

## Interfaces & Ecosystem
- **Abinit**: Primary DFT code
- **Python**: Analysis framework
- **abipy**: Compatible tools


## Advanced Features
- **Allen-Heine-Cardona theory**: Complete temperature dependence
- **Zero-point renormalization**: Quantum effects on band gaps
- **Fan-Migdal self-energy**: Electron-phonon coupling contributions
- **Debye-Waller terms**: Lattice vibration effects
- **Spectral functions**: Full energy-dependent analysis
- **Abinit netCDF**: Direct parsing of Abinit output

## Performance Characteristics
- Post-processing tool: Moderate speed
- Depends on k-point and q-point grids
- Python-based implementation

## Computational Cost
- Abinit DFPT: Dominant cost (external)
- EPC analysis: Minutes to hours
- Scales with system size and grid density
- Overall: DFPT calculations dominate

## Best Practices
- Converge k-point and q-point grids
- Validate against experimental band gap temperature dependence
- Check Fan and Debye-Waller contributions separately
- Use appropriate smearing for spectral functions

## Limitations & Known Constraints
- Abinit-specific
- Requires DFPT calculations
- Python expertise needed
- Limited documentation

## Application Areas
- Semiconductor band gaps
- Temperature-dependent properties
- Superconductivity studies
- Optical properties

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/GkAntonius/ElectronPhononCoupling
2. G. Antonius et al., Phys. Rev. B 92, 085137 (2015)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-3.0)
