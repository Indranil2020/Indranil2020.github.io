# pyFitIt

## Official Resources
- Source Repository: https://github.com/gudasergey/pyFitIt
- Documentation: Included in repository
- License: Open source

## Overview
**pyFitIt** is a Python implementation of the FitIt software for fitting X-ray Absorption Near Edge Structure (XANES) and other spectra. It uses machine learning and automatic component analysis to fit theoretical XANES spectra to experimental data, enabling structural determination from XANES measurements.

**Scientific domain**: XANES fitting, structure determination from spectra  
**Target user community**: Researchers fitting experimental XANES spectra to determine local atomic structure

## Theoretical Methods
- XANES fitting via ML surrogate models
- Automatic component analysis
- Joint convolution fitting
- FEFF-based XANES calculation
- Structural parameter optimization
- Machine learning interpolation
- Grid-based and adaptive sampling

## Capabilities (CRITICAL)
- XANES spectral fitting
- ML-accelerated fitting (surrogate models)
- Automatic component analysis
- Joint convolution fitting
- Structural parameter extraction
- Confidence interval estimation
- Multi-edge fitting
- Grid and adaptive sampling strategies

**Sources**: GitHub repository, J. Synchrotron Rad. 27, 1694 (2020)

## Key Strengths

### ML-Accelerated Fitting:
- Surrogate model replaces expensive FEFF calls
- Orders of magnitude faster fitting
- Enables extensive parameter exploration
- Robust optimization

### Structural Determination:
- Extract atomic positions from XANES
- Determine coordination geometry
- Estimate bond distances
- Confidence intervals on parameters

### Flexible Fitting:
- Multiple fitting strategies
- Custom structural parameters
- Multi-edge constraints
- Experimental resolution convolution

## Inputs & Outputs
- **Input formats**:
  - Experimental XANES spectra
  - FEFF input templates
  - Structural parameter ranges
  
- **Output data types**:
  - Fitted structural parameters
  - Fitted XANES spectra
  - Confidence intervals
  - Residual analysis
  - Component analysis results

## Interfaces & Ecosystem
- **FEFF**: XANES calculation backend
- **Python/Scikit-learn**: ML models
- **NumPy/SciPy**: Optimization
- **Matplotlib**: Visualization

## Performance Characteristics
- **Speed**: Fast with ML surrogate (seconds per fit)
- **Accuracy**: Limited by FEFF accuracy and ML model
- **Parameters**: Typically 3-10 structural parameters
- **Memory**: Low

## Computational Cost
- **ML training**: Hours (FEFF grid calculation)
- **Fitting**: Seconds per spectrum
- **FEFF grid**: Minutes per point × grid size
- **Typical**: Efficient after surrogate training

## Limitations & Known Constraints
- **FEFF-dependent**: Requires FEFF for training data
- **Local structure only**: XANES is local probe
- **Parameter space**: Limited to few structural parameters
- **ML extrapolation**: Poor outside training domain
- **Single-scattering**: Limited multiple scattering accuracy

## Comparison with Other Codes
- **vs XANESNET**: pyFitIt is fitting (inverse), XANESNET is prediction (forward)
- **vs MLXANES**: pyFitIt uses ML for fitting, MLXANES for prediction
- **vs Demeter**: pyFitIt fits structure, Demeter fits EXAFS
- **Unique strength**: ML-accelerated XANES fitting for structural determination, inverse problem solving

## Application Areas

### Catalyst Structure:
- Active site geometry
- Under reaction conditions
- Oxidation state determination
- Adsorbate structure

### Nanoparticle Characterization:
- Size and shape from XANES
- Surface vs bulk structure
- Ligand effects
- Support interactions

### Amorphous Materials:
- Local structure determination
- Bond angle distributions
- Coordination geometry
- Disorder quantification

### Phase Transitions:
- Structural changes across transitions
- Order parameters from XANES
- Time-resolved fitting
- Pressure-dependent structure

## Best Practices

### Grid Construction:
- Cover relevant parameter space
- Use appropriate grid density
- Include edge cases
- Validate ML surrogate quality

### Fitting:
- Start with broad parameter ranges
- Use confidence intervals
- Validate against known structures
- Cross-validate with EXAFS if available

## Community and Support
- Open source on GitHub
- Developed by SMARG group
- Published in J. Synchrotron Rad.
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/gudasergey/pyFitIt
2. A. Martini et al., J. Synchrotron Rad. 27, 1694 (2020)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Published methodology: J. Synchrotron Rad.
- Active development: Maintained
- Specialized strength: ML-accelerated XANES fitting for structural determination
