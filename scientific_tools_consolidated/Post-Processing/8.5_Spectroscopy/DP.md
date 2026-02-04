# DP (Dielectric Properties)

## Official Resources
- Homepage: http://dp-code.org/
- Documentation: http://dp-code.org/documentation/
- Source Repository: Available via website
- License: GNU General Public License
- Developers: ETSF (V. Olevano, L. Reining, G. Siny)

## Overview
DP is a code for calculating the linear response dielectric properties of periodic systems. It uses Time-Dependent Density Functional Theory (TDDFT) in the frequency domain with a plane-wave basis set. It computes the macroscopic dielectric function, Electron Energy Loss Spectra (EELS), and Inelastic X-ray Scattering (IXS) spectra.

**Scientific domain**: Dielectric properties, EELS, optical response
**Target user community**: Materials scientists studying optical and dielectric properties

## Theoretical Methods
- Time-Dependent DFT (TDDFT)
- Random Phase Approximation (RPA)
- Adiabatic LDA (ALDA) kernel
- Long-range corrected (LRC) kernels
- Local field effects

## Capabilities (CRITICAL)
- **Dielectric Function**: Frequency-dependent ε(ω)
- **EELS**: Electron energy loss spectra
- **IXS**: Inelastic X-ray scattering
- **Local Field Effects**: Crystal local fields
- **ABINIT Interface**: Reads WFK files
- **Multiple Kernels**: RPA, ALDA, LRC

**Sources**: DP website, ETSF documentation

## Key Strengths

### TDDFT Implementation:
- Multiple xc kernels
- Local field effects
- Accurate response
- Well-validated

### ABINIT Integration:
- Direct WFK reading
- Consistent workflow
- Plane-wave basis
- ETSF standard

### Open Source:
- GPL licensed
- ETSF developed
- Stable codebase
- Community support

## Inputs & Outputs
- **Input formats**: ABINIT WFK files, DP input file
- **Output data types**: Dielectric function, EELS spectra, IXS spectra

## Performance Characteristics
- Efficient for moderate system sizes
- Scales with k-points and bands
- Parallelized

## Limitations & Known Constraints
- **ABINIT only**: Requires ABINIT wavefunctions
- **Legacy code**: Stable but less active development
- **Documentation**: Could be more extensive
- **Learning curve**: TDDFT concepts required

## Comparison with Other Tools
- **vs Yambo**: DP simpler, Yambo more features
- **vs GPAW**: Different implementations
- **vs exciting**: DP plane-wave, exciting LAPW
- **Unique strength**: ETSF standard, ABINIT integration

## Application Areas
- Optical properties of solids
- Plasmon analysis
- EELS simulation
- Dielectric screening

## Best Practices
- Converge k-points and bands
- Test local field effects
- Compare RPA vs ALDA
- Validate with experiment

## Community and Support
- ETSF development
- GPL licensed
- Stable/legacy status
- Academic support

## Verification & Sources
**Primary sources**:
1. Homepage: http://dp-code.org/
2. ETSF documentation

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACCESSIBLE
- Source: OPEN (GPL)
- Status: Stable/Legacy
- Method: TDDFT for dielectric properties
