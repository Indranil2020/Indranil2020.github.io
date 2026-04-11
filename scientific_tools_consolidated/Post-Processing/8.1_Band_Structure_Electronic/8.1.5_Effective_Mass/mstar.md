# mstar

## Official Resources
- Source Repository: https://github.com/rubel75/mstar
- Documentation: Included in repository
- License: Open source

## Overview
**mstar** is a tool for calculating effective masses with DFT using perturbation theory. It computes conductivity, density-of-states, and cyclotron effective masses from WIEN2k band structures using the k·p perturbation theory approach, providing more accurate effective masses than finite-difference methods.

**Scientific domain**: Effective mass calculation, semiconductor band structure analysis  
**Target user community**: Researchers studying carrier effective masses in semiconductors and insulators

## Theoretical Methods
- k·p perturbation theory for effective masses
- Conductivity effective mass (m*_c)
- Density-of-states effective mass (m*_dos)
- Cyclotron effective mass (m*_cyc)
- Principal components of inverse mass tensor
- WIEN2k band structure input

## Capabilities (CRITICAL)
- Conductivity effective mass calculation
- Density-of-states effective mass
- Cyclotron effective mass
- Inverse mass tensor principal components
- Perturbation theory approach (more accurate than finite differences)
- WIEN2k interface

**Sources**: GitHub repository, Comp. Phys. Commun.

## Key Strengths

### Perturbation Theory:
- More accurate than finite-difference methods
- Systematic convergence
- No numerical differentiation errors
- Well-defined at band extrema

### Multiple Mass Definitions:
- Conductivity mass (transport)
- DOS mass (thermodynamics)
- Cyclotron mass (cyclotron resonance)
- Anisotropic mass tensor

### WIEN2k Integration:
- Direct interface with WIEN2k
- Uses same basis sets
- Consistent calculation flow

## Inputs & Outputs
- **Input formats**:
  - WIEN2k band structure output
  - k-point specifications
  
- **Output data types**:
  - Conductivity effective mass (m0/m*_c)
  - DOS effective mass
  - Cyclotron effective mass
  - Principal components of inverse mass tensor

## Interfaces & Ecosystem
- **WIEN2k**: DFT backend
- **Fortran**: Core computation

## Performance Characteristics
- **Speed**: Fast (post-processing)
- **Accuracy**: High (perturbation theory)
- **System size**: Limited by WIEN2k
- **Memory**: Low

## Computational Cost
- **Effective mass**: Seconds
- **WIEN2k pre-requisite**: Hours (separate)
- **Typical**: Very efficient

## Limitations & Known Constraints
- **WIEN2k only**: No VASP or QE support
- **Band extrema only**: Requires identified extrema
- **Perturbation theory**: Limited to parabolic regions
- **Documentation**: Limited

## Comparison with Other Codes
- **vs effmass**: mstar uses perturbation theory, effmass uses finite differences from VASP
- **vs emc**: mstar is WIEN2k, emc is VASP/QE finite differences
- **vs Effective-mass-fitting**: mstar is perturbation theory, fitting is polynomial
- **Unique strength**: Perturbation theory effective masses from WIEN2k, multiple mass definitions

## Application Areas

### Semiconductor Physics:
- Carrier effective mass for transport
- DOS mass for thermodynamic properties
- Cyclotron mass for magneto-optics
- Anisotropic mass for device modeling

### Thermoelectric Materials:
- DOS mass for Seebeck coefficient
- Conductivity mass for electrical conductivity
- Mass anisotropy for direction-dependent transport
- Effective mass optimization

### Optoelectronic Materials:
- Reduced mass for exciton binding
- Effective mass for optical absorption
- Carrier mobility estimation
- Band structure engineering

## Best Practices

### WIEN2k Setup:
- Use well-converged SCF calculation
- Adequate k-point density near extrema
- Include spin-orbit coupling if needed
- Use consistent settings

### Mass Calculation:
- Identify band extrema correctly
- Use sufficient k-points near extrema
- Compare perturbation vs finite difference
- Validate against experimental cyclotron resonance

## Community and Support
- Open source on GitHub
- Developed by O. Rubel
- Published in Comp. Phys. Commun.
- Research code

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/rubel75/mstar
2. O. Rubel, F. Tran, X. Rocquefelte, and P. Blaha, Comp. Phys. Commun. (related)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Published methodology: Comp. Phys. Commun.
- Specialized strength: Perturbation theory effective masses from WIEN2k, multiple mass definitions
