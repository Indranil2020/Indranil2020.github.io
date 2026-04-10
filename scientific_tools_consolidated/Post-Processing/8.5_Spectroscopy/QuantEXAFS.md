# QuantEXAFS

## Official Resources
- Source Repository: https://github.com/kul-group/QuantEXAFS
- Documentation: Included in repository
- License: Open source

## Overview
**QuantEXAFS** is a Python-based toolkit for automated fitting of Extended X-ray Absorption Fine Structure (EXAFS) data using X-ray Larch modules. It combines DFT-optimized structure databases with automated EXAFS fitting workflows, enabling quantitative structural analysis from EXAFS measurements.

**Scientific domain**: EXAFS analysis, X-ray absorption spectroscopy, structural fitting  
**Target user community**: Researchers performing quantitative EXAFS analysis with DFT structural models

## Theoretical Methods
- EXAFS theory (scattering formalism)
- FEFF scattering path calculation
- DFT structure optimization
- Nonlinear least-squares fitting
- Artemis/Larch fitting engine
- Multiple scattering paths

## Capabilities (CRITICAL)
- Automated EXAFS fitting
- DFT structure database integration
- Multiple scattering path analysis
- Fitting of coordination numbers
- Fitting of Debye-Waller factors
- Fitting of bond distances
- Batch fitting of multiple spectra
- ASE database format support

**Sources**: GitHub repository

## Key Strengths

### Automated Workflow:
- No manual fitting steps
- Reproducible results
- Batch processing
- Systematic parameter exploration

### DFT Integration:
- DFT-optimized structures as starting models
- ASE database format
- Consistent structure-spectra workflow
- High-quality structural models

### Larch Integration:
- Uses mature XAS analysis library
- Well-tested fitting algorithms
- Standard EXAFS methodology
- Community-validated

## Inputs & Outputs
- **Input formats**:
  - EXAFS data files
  - ASE database of DFT structures
  - FEFF calculation results
  - Fitting parameter files
  
- **Output data types**:
  - Fitted structural parameters
  - Fitted EXAFS spectra
  - Residuals and R-factors
  - Coordination numbers
  - Bond distances and disorder

## Interfaces & Ecosystem
- **Larch**: XAS analysis library
- **ASE**: Structure database management
- **FEFF**: Scattering path calculation
- **DFT codes**: Structure optimization (VASP, QE, etc.)

## Performance Characteristics
- **Speed**: Fast (fitting is seconds, FEFF is minutes)
- **Accuracy**: Depends on model quality
- **System size**: Any size (EXAFS is local probe)
- **Memory**: Low

## Computational Cost
- **Fitting**: Seconds per spectrum
- **FEFF paths**: Minutes per absorber
- **DFT structures**: Hours (pre-requisite)
- **Typical**: Very efficient fitting step

## Limitations & Known Constraints
- **EXAFS only**: No XANES fitting
- **Larch dependency**: Requires Larch installation
- **DFT pre-requisite**: Need DFT-optimized structures
- **Documentation**: Limited
- **Local structure only**: EXAFS probes local environment

## Comparison with Other Codes
- **vs Demeter/Athena**: QuantEXAFS adds DFT database integration
- **vs Larch**: QuantEXAFS automates fitting workflow
- **vs FEFF**: QuantEXAFS uses FEFF for paths, adds fitting
- **Unique strength**: Automated EXAFS fitting with DFT structure database integration

## Application Areas

### Nanoparticle Structure:
- Size-dependent structure
- Surface vs bulk coordination
- Shape determination
- Ligand effects

### Amorphous Materials:
- Short-range order
- Bond distance distributions
- Coordination number analysis
- Structural modeling

### Catalysis:
- Active site structure
- Under operating conditions
- Catalyst degradation
- Support effects

### Battery Materials:
- Local structure changes
- Phase transitions
- Cation disorder
- Redox processes

## Best Practices

### DFT Structures:
- Use well-converged geometries
- Include relevant structural models
- Consider disorder and defects
- Validate against known structures

### EXAFS Fitting:
- Use appropriate k-range
- Include sufficient R-range
- Test fitting stability
- Report uncertainties

## Community and Support
- Open source on GitHub
- Developed at KU Leuven
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/kul-group/QuantEXAFS
2. Related publications from KU Leuven

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Research code
- Specialized strength: Automated EXAFS fitting with DFT structure database integration
