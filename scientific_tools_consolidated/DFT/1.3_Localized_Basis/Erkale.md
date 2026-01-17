# Erkale

## Official Resources
- Homepage: https://github.com/susilehtola/erkale
- Documentation: https://github.com/susilehtola/erkale/wiki
- Source Repository: https://github.com/susilehtola/erkale
- License: GNU General Public License v2.0

## Overview
Erkale is an open-source quantum chemistry program for Hartree-Fock and density functional theory calculations using Gaussian basis sets. Originally developed at the University of Helsinki, it focuses on computing X-ray properties including ground-state electron momentum densities, Compton profiles, X-ray absorption spectra (XAS), and X-ray Raman scattering spectra. Erkale has evolved to include advanced capabilities in basis set development and self-interaction corrected DFT.

**Scientific domain**: Molecules, X-ray spectroscopy, core-level physics, basis set development  
**Target user community**: Researchers studying X-ray properties, spectroscopy, and those needing SIC-DFT or basis set optimization tools

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Density Functional Theory (DFT)
- Gaussian Type Orbitals (GTOs) of arbitrary angular momentum
- Exchange-correlation via LibXC (600+ functionals)
- Self-Interaction Correction (SIC-DFT)
- Orbital localization methods (Foster-Boys, Pipek-Mezey, etc.)
- Core-level spectroscopy methods
- Basis set optimization algorithms

## Capabilities (CRITICAL)
- Ground-state electronic structure
- X-ray absorption spectroscopy (XAS) simulation
- X-ray Raman scattering spectra
- Compton profiles
- Electron momentum densities
- Core-level excitations (K-edge, L-edge)
- Self-interaction corrected DFT
- Basis set completeness optimization
- Multiple orbital localization schemes
- Transition potential method
- Full core-hole approximation

**Sources**: GitHub repository, University of Helsinki, published papers

## Key Strengths

### X-ray Spectroscopy Focus:
- Native XAS simulation capabilities
- X-ray Raman scattering
- Compton profile calculations
- Core-level physics expertise
- Direct comparison with synchrotron experiments

### Self-Interaction Correction:
- Perdew-Zunger SIC implementation
- Improved orbital energies
- Better description of localized states
- Corrected band gaps

### Basis Set Development:
- Automatic basis set optimization
- Completeness-optimized basis sets
- Angular momentum extensions
- Contraction optimization

### Modern Implementation:
- Object-oriented C++ design
- LibXC integration (600+ functionals)
- ADIIS/Broyden convergence accelerators
- Easy to understand and extend

## Inputs & Outputs
- **Input formats**:
  - Native Erkale input files
  - XYZ coordinates
  - Basis set specifications (Gaussian format)
  
- **Output data types**:
  - Total energies
  - Orbital energies and coefficients
  - XAS spectra
  - Compton profiles
  - Electron momentum densities
  - Localized orbitals

## Interfaces & Ecosystem
- **LibXC integration**:
  - Access to 600+ density functionals
  - LDA, GGA, meta-GGA, hybrid functionals
  - Range-separated hybrids
  
- **Basis set libraries**:
  - Standard Gaussian basis formats
  - Basis Set Exchange compatibility
  - Custom basis optimization

- **Visualization**:
  - Molden format output
  - Cube file generation
  - Standard plotting tools

## Advanced Features

### Core-Level Spectroscopy:
- Full core-hole approximation
- Transition potential method
- Core-valence separation
- Element-specific probing
- Comparison with experimental spectra

### Orbital Localization:
- Foster-Boys localization
- Pipek-Mezey localization
- Edmiston-Ruedenberg
- Fourth-moment methods
- Intrinsic atomic orbitals

### Basis Set Optimization:
- Completeness profiles
- Exponent optimization
- Polarization function addition
- Contraction schemes
- Element-specific tuning

### SIC-DFT:
- Perdew-Zunger formulation
- Improved ionization potentials
- Better charge localization
- Reduced self-interaction error

## Performance Characteristics
- **Speed**: Efficient for medium-sized systems
- **Accuracy**: High accuracy for spectroscopy
- **System size**: Molecules up to ~100 atoms
- **Memory**: Standard Gaussian code requirements
- **Parallelization**: OpenMP threading

## Computational Cost
- **DFT/HF**: Standard Gaussian scaling
- **XAS**: Additional cost for core-hole calculations
- **SIC**: 2-3x overhead per iteration
- **Typical**: Desktop calculations for most molecules
- **Large basis**: Feasible with thousands of functions

## Limitations & Known Constraints
- **Periodicity**: Molecular only (no periodic systems)
- **System size**: Best for small to medium molecules
- **Forces**: Limited geometry optimization
- **User base**: Specialized (X-ray spectroscopy)
- **Documentation**: Academic-level
- **Dynamics**: No molecular dynamics

## Comparison with Other Codes
- **vs Gaussian/ORCA**: Erkale specialized for X-ray, general codes broader
- **vs FLOSIC**: Both have SIC-DFT, different implementations
- **vs StoBe**: Both X-ray focused, different methodologies
- **Unique strength**: X-ray spectroscopy, SIC-DFT, basis set development, open-source

## Application Areas

### X-ray Spectroscopy:
- K-edge XANES simulation
- L-edge spectra
- X-ray Raman scattering
- Comparison with synchrotron data
- Element-specific probing

### Core-Level Physics:
- Core ionization potentials
- Chemical shifts
- Core-hole effects
- Auger processes

### Electronic Structure:
- Ground-state calculations
- Orbital localization analysis
- Electron density analysis
- Bonding characterization

### Basis Set Research:
- Completeness optimization
- New basis set development
- Polarization function design
- Contraction schemes

## Best Practices

### XAS Calculations:
- Use appropriate core-hole approximation
- Converge basis set for core region
- Compare with experimental calibration
- Account for relativistic effects for heavy elements

### Basis Set Selection:
- Start with standard basis (cc-pVTZ)
- Add diffuse functions for Rydberg states
- Test core-region completeness
- Document basis choice

### SIC Calculations:
- Start from converged standard DFT
- Monitor SIC energy convergence
- Check orbital localization
- Compare with experimental IPs

### SCF Convergence:
- Use ADIIS for difficult cases
- Monitor energy convergence
- Level shifting if needed

## Community and Support
- Open source GPL v2
- GitHub repository
- Academic publications
- Author-maintained (S. Lehtola)
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/susilehtola/erkale
2. S. Lehtola et al., J. Comput. Chem. 33, 1572 (2012)
3. S. Lehtola, J. Chem. Theory Comput. publications

**Secondary sources**:
1. X-ray spectroscopy literature
2. Basis set optimization papers
3. SIC-DFT methodology papers

**Confidence**: VERIFIED - Active GitHub, published methodology

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL v2)
- Academic use: Published applications
- Documentation: Wiki and papers
- Active development: Recent commits
- Specialty: X-ray spectroscopy, SIC-DFT, basis set development
