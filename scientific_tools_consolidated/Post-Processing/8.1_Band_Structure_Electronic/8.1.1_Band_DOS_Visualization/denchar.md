# denchar

## Official Resources
- **SIESTA**: https://siesta-project.org/
- **Part of**: SIESTA utilities
- **License**: GPL

## Overview
denchar is a SIESTA utility for charge density and wavefunction plotting, providing 2D/3D visualization of electron density and wavefunctions from SIESTA calculations.

**Scientific domain**: SIESTA post-processing, charge density visualization
**Target user community**: SIESTA users

## Capabilities (CRITICAL)
- **Charge Density**: 2D/3D electron density plots
- **Wavefunctions**: Orbital visualization
- **Cube Files**: Output in Gaussian cube format
- **Slices**: 2D cuts through 3D data

## Key Strengths
- SIESTA native utility
- 2D and 3D visualization
- Cube file output
- Part of SIESTA distribution

## Inputs & Outputs
- **Input formats**: SIESTA output files (.RHO, .WFS)
- **Output data types**: Cube files, 2D data files

## Installation
Included with SIESTA distribution. Build from SIESTA utilities:
```bash
cd siesta/Util/Denchar
make
```

## Limitations & Known Constraints
- **SIESTA-specific**: Only works with SIESTA output
- **Compilation**: Requires building with SIESTA
- **Output format**: Limited to cube files and 2D slices

## Comparison with Other Tools
- **vs VESTA**: denchar SIESTA-native, VESTA general visualization
- **vs c2x**: Both can produce cube files, different ecosystems
- **Unique strength**: Native SIESTA charge/wavefunction visualization

## Verification & Sources
**Confidence**: VERIFIED - Part of SIESTA

**Verification status**: ✅ VERIFIED
- Part of SIESTA distribution
- Target Code: SIESTA
