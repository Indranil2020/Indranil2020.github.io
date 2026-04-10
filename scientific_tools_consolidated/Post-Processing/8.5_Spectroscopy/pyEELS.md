# pyEELS

## Official Resources
- Source Repository: https://github.com/sindrebilden/pyeels
- Documentation: Included in repository
- License: Open source

## Overview
**pyEELS** is a Python simulation package for constructing Electron Energy Loss Spectroscopy (EELS) spectra from model band structures. It simulates EELS based on model band structures created using PythTB (Python Tight Binding) or parabolic band models, enabling rapid spectral prediction for materials characterization.

**Scientific domain**: Electron energy loss spectroscopy, band structure-based simulation  
**Target user community**: Researchers simulating EELS spectra for comparison with STEM-EELS experiments

## Theoretical Methods
- Electron energy loss spectroscopy (EELS) theory
- Dielectric function formalism
- Model band structures (PythTB, parabolic)
- Single and double differential cross-sections
- Momentum-resolved EELS
- Core-loss and low-loss EELS

## Capabilities (CRITICAL)
- Low-loss EELS simulation
- Core-loss EELS simulation
- Model band structure input
- PythTB integration
- Parabolic band model
- Momentum-resolved spectra
- Thickness-dependent spectra
- Comparison with experimental EELS

**Sources**: GitHub repository, University of Oslo thesis

## Key Strengths

### Model-Based:
- Fast simulation from model bands
- No DFT calculation required
- Rapid exploration of parameter space
- Educational tool for EELS understanding

### PythTB Integration:
- Tight-binding band structures
- Customizable Hamiltonians
- Well-established TB package
- Flexible model construction

### Python Implementation:
- Easy to use and modify
- Jupyter notebook compatible
- Visualization included
- Extensible framework

## Inputs & Outputs
- **Input formats**:
  - PythTB model definitions
  - Parabolic band parameters
  - EELS geometry parameters
  
- **Output data types**:
  - EELS spectra (energy loss vs intensity)
  - Momentum-resolved EELS
  - Dielectric function
  - Cross-sections

## Interfaces & Ecosystem
- **PythTB**: Tight-binding band structure
- **NumPy**: Numerical computation
- **Matplotlib**: Visualization
- **Python**: Scripting

## Performance Characteristics
- **Speed**: Very fast (model-based)
- **Accuracy**: Qualitative, depends on model quality
- **System size**: Any (model-based)
- **Memory**: Low

## Computational Cost
- **Per spectrum**: Seconds
- **Parameter scan**: Minutes
- **Typical**: Very efficient

## Limitations & Known Constraints
- **Model-based**: Not ab initio
- **Qualitative**: Not quantitative accuracy
- **No core-hole effects**: Missing for core-loss
- **No multiplet effects**: Single-particle
- **Limited documentation**: Research code

## Comparison with Other Codes
- **vs FEFF**: pyEELS is model-based, FEFF is ab initio
- **vs HyperSpy**: pyEELS simulates, HyperSpy analyzes experimental data
- **vs turboEELS**: pyEELS is model-based, turboEELS is TDDFT-based
- **Unique strength**: Fast EELS simulation from model band structures, PythTB integration

## Application Areas

### Low-Loss EELS:
- Plasmon spectroscopy
- Band gap determination
- Dielectric function
- Interband transitions

### Core-Loss EELS:
- Edge onset modeling
- White line intensity
- Chemical shift estimation
- Fingerprinting

### STEM-EELS:
- Spectrum image simulation
- Line scan simulation
- Thickness dependence
- Comparison with experiment

## Best Practices

### Model Selection:
- Choose appropriate band structure model
- Validate against DFT if available
- Test sensitivity to model parameters
- Compare with experimental EELS

### Spectral Analysis:
- Account for experimental resolution
- Include zero-loss peak for low-loss
- Consider thickness effects
- Use appropriate broadening

## Community and Support
- Open source on GitHub
- Developed at University of Oslo
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/sindrebilden/pyeels
2. S. Bilden, Master thesis, University of Oslo

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Limited
- Active development: Research code
- Specialized strength: Fast EELS simulation from model band structures, PythTB integration
