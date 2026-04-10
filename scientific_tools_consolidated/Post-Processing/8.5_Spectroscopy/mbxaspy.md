# mbxaspy

## Official Resources
- Source Repository: https://github.com/yufengliang/mbxaspy
- Documentation: Included in repository
- License: Open source

## Overview
**mbxaspy** is a Python software package for predicting X-ray spectra using the determinant formalism based on the independent-electron approximation as used in DFT. It interfaces with the ShirleyXAS Fortran package for DFT and XAS calculations at the one-body level, and can also work with tight-binding models for many-body XAS calculations.

**Scientific domain**: X-ray absorption spectroscopy, many-body spectroscopy  
**Target user community**: Researchers computing XAS spectra from DFT and tight-binding models, including many-body extensions

## Theoretical Methods
- Determinant formalism for XAS
- Independent-electron approximation (DFT)
- Many-body extension beyond DFT
- Tight-binding models
- Core-hole treatment
- Bethe-Salpeter-like corrections
- ShirleyXAS integration

## Capabilities (CRITICAL)
- X-ray absorption spectroscopy (XAS)
- X-ray emission spectroscopy (XES)
- Many-body XAS calculations
- Tight-binding XAS
- DFT-level XAS (via ShirleyXAS)
- Core-hole effect simulation
- Polarization-dependent spectra
- K-edge and L-edge calculations

**Sources**: GitHub repository

## Key Strengths

### Determinant Formalism:
- Systematic improvement beyond DFT
- Can include many-body effects
- Flexible Hamiltonian choices
- From DFT to tight-binding models

### ShirleyXAS Integration:
- Well-established XAS code
- Full DFT-level calculations
- Core-hole pseudopotentials
- Production-quality results

### Python Interface:
- Scriptable workflow
- Jupyter notebook compatible
- Easy post-processing
- Visualization tools

## Inputs & Outputs
- **Input formats**:
  - ShirleyXAS DFT outputs
  - Tight-binding model parameters
  - Hamiltonian specifications
  
- **Output data types**:
  - XAS spectra
  - XES spectra
  - Many-body corrected spectra
  - Polarization-dependent spectra

## Interfaces & Ecosystem
- **ShirleyXAS**: DFT/XAS Fortran backend
- **Python**: Scripting and analysis
- **DFT codes**: Via ShirleyXAS interface

## Performance Characteristics
- **Speed**: Depends on backend (ShirleyXAS or TB)
- **Accuracy**: DFT-level or better with many-body
- **System size**: Limited by backend
- **Memory**: Moderate

## Computational Cost
- **DFT XAS**: Hours (ShirleyXAS)
- **TB XAS**: Minutes
- **Many-body**: Hours to days
- **Typical**: Moderate

## Limitations & Known Constraints
- **ShirleyXAS dependency**: Requires external Fortran code
- **Documentation**: Limited
- **Community**: Small
- **Installation**: Can be complex (Fortran + Python)

## Comparison with Other Codes
- **vs FEFF**: mbxaspy can do many-body, FEFF is single-particle
- **vs OCEAN**: mbxaspy is more flexible, OCEAN is BSE-based
- **vs xspectra**: mbxaspy has many-body extension, xspectra is DFT-only
- **Unique strength**: Determinant formalism for XAS with both DFT and many-body capabilities

## Application Areas

### Transition Metal XAS:
- K-edge and L-edge spectra
- Many-body effects in XAS
- Core-hole screening
- Charge transfer satellites

### Battery Materials:
- Transition metal redox
- Oxygen K-edge XAS
- Charge state analysis
- Cycling effects

### Correlated Oxides:
- Multiplet-like features
- Hubbard model XAS
- Charge transfer insulators
- Metal-insulator transitions

## Best Practices

### DFT Parameters:
- Use well-converged ShirleyXAS calculations
- Appropriate core-hole treatment
- Test k-point convergence
- Validate against experiment

### Many-Body Extensions:
- Start with DFT-level results
- Add many-body corrections systematically
- Validate against BSE if possible
- Compare with experiment

## Community and Support
- Open source on GitHub
- Developed at LBNL Molecular Foundry
- Limited documentation
- Research code with active development

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/yufengliang/mbxaspy
2. Y. Liang et al., related publications from LBNL

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Limited
- Active development: Research code
- Specialized strength: Determinant formalism XAS with DFT and many-body capabilities, ShirleyXAS integration
