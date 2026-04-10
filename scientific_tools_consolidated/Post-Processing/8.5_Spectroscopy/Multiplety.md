# Multiplety

## Official Resources
- Source Repository: https://github.com/gfabbris/multiplety
- License: Open source

## Overview
**Multiplety** is a Python package for multiplet calculations of X-ray absorption (XAS) and resonant inelastic X-ray scattering (RIXS) spectra using the Cowan's atomic code and Racer programs. It provides a Jupyter notebook-based interface for setting up and running multiplet calculations for transition metal and rare-earth systems.

**Scientific domain**: Multiplet X-ray spectroscopy, atomic physics  
**Target user community**: Researchers studying core-level spectra of correlated transition metal and rare-earth systems

## Theoretical Methods
- Cowan's atomic code (Hartree-Fock with relativistic corrections)
- Racer code for transition matrices
- Crystal field theory
- Configuration interaction (CI)
- Spin-orbit coupling
- Slater-Condon parameter scaling
- Multiplet theory

## Capabilities (CRITICAL)
- X-ray absorption spectroscopy (XAS) multiplet calculations
- RIXS multiplet calculations
- XAS dichroism (XMCD, XMLD)
- Temperature-dependent spectra
- Polarization-dependent spectra
- Interactive Jupyter notebook interface
- Automatic parameter setup from Cowan's code
- Custom crystal field splitting
- Slater-Condon parameter scaling

**Sources**: GitHub repository, Cowan's code documentation

## Key Strengths

### Cowan's Code Integration:
- Well-established atomic multiplet code
- Accurate radial integrals
- Relativistic corrections included
- Systematic parameter scaling
- Decades of validation

### User-Friendly Interface:
- Jupyter notebook workflow
- Interactive parameter adjustment
- Visual output
- Step-by-step tutorials
- Python scripting

### Spectroscopy Coverage:
- XAS with multiplet structure
- RIXS with full multiplet treatment
- Dichroism calculations
- Temperature effects
- Polarization dependence

## Inputs & Outputs
- **Input formats**:
  - Python/Jupyter notebooks
  - Atomic configuration specifications
  - Crystal field parameters
  - Slater-Condon scaling factors
  
- **Output data types**:
  - XAS spectra
  - RIXS maps
  - Dichroism spectra
  - Energy level diagrams
  - Transition matrices

## Interfaces & Ecosystem
- **Cowan's code**: Required external dependency
- **Racer**: Required for transition matrices
- **Jupyter**: Interactive interface
- **Matplotlib**: Visualization

## Performance Characteristics
- **Speed**: Fast for single-atom multiplet calculations
- **Accuracy**: Good for atomic multiplet structure
- **System size**: Single atom/impurity models
- **Memory**: Low

## Computational Cost
- **XAS**: Seconds to minutes
- **RIXS**: Minutes
- **Typical**: Very fast for atomic calculations

## Limitations & Known Constraints
- **Requires Cowan's code**: External dependency needed
- **Single atom**: No ligand/cluster effects natively
- **No DFT integration**: Standalone multiplet only
- **Installation**: Cowan's code setup can be complex
- **Documentation**: Limited

## Comparison with Other Codes
- **vs Quanty**: Multiplety uses Cowan's code, Quanty is Lua-based
- **vs CTM4XAS**: Multiplety is Python, CTM4XAS is MATLAB
- **vs Crispy**: Multiplety is standalone, Crispy wraps Quanty
- **Unique strength**: Direct Cowan's code integration, Python/Jupyter interface, simple setup

## Application Areas

### Transition Metal L-edges:
- 2p→3d XAS multiplet structure
- Crystal field splitting analysis
- Oxidation state determination
- Spin-state characterization

### Rare-Earth M-edges:
- 3d→4f XAS multiplet structure
- Hund's rule coupling
- Crystal field effects
- Valence determination

### RIXS of Correlated Systems:
- d-d excitations
- Charge transfer excitations
- Parametric studies
- Temperature dependence

## Best Practices

### Parameter Calibration:
- Scale Slater-Condon parameters (typically 70-85%)
- Validate against experimental spectra
- Use consistent scaling across edges
- Compare with charge transfer models

### Crystal Field Setup:
- Use appropriate symmetry
- Start with DFT crystal field parameters
- Validate splitting patterns
- Consider low-symmetry distortions

## Community and Support
- Open source on GitHub
- Limited but active development
- Jupyter notebook examples provided
- Based on well-established Cowan's code

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/gfabbris/multiplety
2. Cowan's code: R. D. Cowan, "The Theory of Atomic Structure and Spectra" (1981)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Jupyter notebooks
- Active development: Maintained
- Specialized strength: Cowan's code integration for multiplet XAS/RIXS, Python/Jupyter interface
