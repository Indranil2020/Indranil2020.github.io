# latticeDynamics

## Official Resources
- Homepage: https://github.com/jgwillingham/latticeDynamics
- Source Repository: https://github.com/jgwillingham/latticeDynamics
- License: MIT License

## Overview
latticeDynamics is a Python package for lattice dynamics calculations using rigid ion and shell models. It provides tools for constructing dynamical matrices, calculating phonon dispersions, and analyzing vibrational properties of crystals.

**Scientific domain**: Lattice dynamics, rigid ion models, phonon calculations  
**Target user community**: Researchers studying phonon properties with empirical potentials

## Theoretical Methods
- Rigid ion model
- Shell model (in development)
- Coulomb interactions (Ewald summation)
- Short-range potentials
- Dynamical matrix construction
- Phonon dispersion relations

## Capabilities (CRITICAL)
- Lattice and Slab structure classes
- Rigid ion model implementation
- Coulomb interaction handling
- Dynamical matrix diagonalization
- Phonon dispersion calculation
- Surface phonon calculations
- Symmetry classification (planned)

## Key Strengths

### Model Flexibility:
- Rigid ion models
- Coulomb interactions
- Short-range potentials
- Extensible framework

### Educational Value:
- Clear Python implementation
- Well-documented code
- Good for learning
- Modifiable

### Surface Support:
- Slab calculations
- Surface phonons
- Interface studies

## Inputs & Outputs
- **Input formats**:
  - Crystal structure definitions
  - Potential parameters
  - Q-point specifications
  
- **Output data types**:
  - Phonon frequencies
  - Eigenvectors
  - Dispersion curves
  - DOS

## Interfaces & Ecosystem
- **Crystals**: Structure handling
- **NumPy**: Numerical operations
- **Matplotlib**: Plotting
- **Python**: Pure Python implementation

## Limitations & Known Constraints
- Empirical potential focus
- Limited to supported models
- Smaller user base
- Some features in development

## Comparison with Other Codes
- **vs Phonopy**: latticeDynamics uses empirical potentials; Phonopy uses DFT forces
- **vs OpenPhonon**: Similar scope; latticeDynamics is Python-native
- **Unique strength**: Clean Python implementation for rigid ion/shell models

## Best Practices

### Model Setup:
- Use appropriate potential parameters
- Validate against experimental data
- Check Coulomb convergence
- Test different cutoffs

### Calculations:
- Use sufficient q-point sampling
- Check acoustic sum rules
- Validate with known results

## Application Areas
- Ionic crystals
- Oxide materials
- Surface phonons
- Educational purposes
- Model development

## Community and Support
- Open-source MIT License
- GitHub repository
- Educational focus
- Python-native

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/jgwillingham/latticeDynamics

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Educational tool
