# kspy-tddft

## Official Resources
- Homepage: https://github.com/pwborthwick/kspy-tddft
- Source Repository: https://github.com/pwborthwick/kspy-tddft
- License: Open Source (Educational)

## Overview
kspy-tddft is an educational Python implementation of Time-Dependent Density Functional Theory (TDDFT) in the Tamm-Dancoff Approximation (TDA). It provides both Linear-Response TDDFT (LR-TDDFT) for computing singlet and triplet excited states, and Real-Time TDDFT (RT-TDDFT) for time propagation simulations. The code is designed for teaching and understanding the fundamental machinery of TDDFT as used in production quantum chemistry codes.

**Scientific domain**: Molecular excited states, UV-Vis spectroscopy, response properties  
**Target user community**: Students and researchers learning TDDFT methodology, educators in computational chemistry

## Theoretical Methods
- Density Functional Theory (DFT)
- Time-Dependent DFT in Tamm-Dancoff Approximation (TDA)
- Linear-Response TDDFT (LR-TDDFT)
- Real-Time TDDFT (RT-TDDFT)
- Slater LDA exchange functional
- RPA correlation functional
- Spin-polarized functionals
- Casida equation formalism

## Capabilities
- Ground-state DFT calculations
- Singlet excited state energies
- Triplet excited state energies
- Electric transition dipoles (length and velocity gauges)
- Magnetic transition dipoles (length gauge)
- Oscillator strengths
- Rotary strengths (both gauges)
- Transition natural orbital (TNO) analysis
- Absorption spectrum plotting
- Real-time electron dynamics

## Key Strengths

### Educational Design:
- Clear, readable Python implementation
- Step-by-step documentation
- Direct comparison with theory
- Minimal dependencies
- Well-documented algorithms

### Complete LR-TDDFT Implementation:
- Singlet and triplet states
- Multiple response properties
- Both length and velocity gauges
- Transition natural orbitals
- Spectrum visualization

### Dual TDDFT Approaches:
- Linear-response for excitation energies
- Real-time for dynamics
- Same code base for comparison

## Inputs & Outputs
- **Input formats**:
  - Hard-coded molecular geometry (H2O default)
  - STO-3G basis set (psi4 format from BSE)
  - Python configuration
  
- **Output data types**:
  - Excitation energies
  - Transition dipoles
  - Oscillator strengths
  - Rotary strengths
  - Transition natural orbitals
  - Plotted spectra

## Interfaces & Ecosystem
- **Dependencies**:
  - NumPy for linear algebra
  - Matplotlib for plotting
  - Standard Python libraries

- **Related projects**:
  - Part of kspy educational quantum chemistry suite
  - Related to harpy HF/post-HF implementation

## Theoretical Background
The code solves the Casida equation for excitation energies, following the formalism described in:
- M.E. Casida, J. Mol. Struct.: THEOCHEM 914 (2009) 3–18
- A. Dreuw and M. Head-Gordon, Chem. Rev. 2005, 105, 4009−4037

The coupling matrix requires second derivatives of the exchange-correlation energy, implemented analytically for Slater LDA exchange and RPA correlation functionals.

## Performance Characteristics
- **Speed**: Educational focus, not optimized for production
- **System size**: Small molecules (H2O, similar)
- **Memory**: Minimal requirements
- **Purpose**: Teaching and prototyping

## Limitations & Known Constraints
- **Basis sets**: Limited to STO-3G (hard-coded)
- **Molecules**: H2O hard-coded as default
- **Functionals**: Only Slater LDA + RPA
- **Performance**: Not production-optimized
- **Features**: No gradients or geometry optimization

## Comparison with Other Codes
- **vs PySCF**: kspy-tddft is purely educational, PySCF is production
- **vs Gaussian/ORCA**: kspy-tddft shows inner workings, production codes are black boxes
- **vs Working_TDDFT**: Both educational, kspy-tddft is functional code
- **Unique strength**: Complete, working LR-TDDFT for learning purposes

## Application Areas
- Teaching TDDFT methodology
- Understanding excited state calculations
- Learning response property theory
- Prototyping new TDDFT methods
- Validating theoretical derivations

## Best Practices
- Use for understanding, not production
- Study code alongside Casida/Dreuw-Head-Gordon papers
- Compare results with production codes
- Extend to new functionals for learning

## Community and Support
- Open-source on GitHub
- Part of educational kspy suite
- Python 100%
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/pwborthwick/kspy-tddft
2. Associated kspy quantum chemistry suite

**Confidence**: VERIFIED - Active GitHub repository with documentation

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Documentation: README and results.md
- Purpose: Educational
- Language: Python 100%
