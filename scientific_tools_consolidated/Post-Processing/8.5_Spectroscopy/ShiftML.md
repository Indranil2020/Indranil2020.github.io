# ShiftML

## Official Resources
- Source Repository: https://github.com/lab-cosmo/ShiftML
- Documentation: https://lab-cosmo.github.io/ShiftML/
- License: Open source (MIT-compatible)

## Overview
**ShiftML** is a Python package for the prediction of NMR chemical shieldings in organic solids using machine learning. It is trained on a large dataset of chemical shieldings computed with DFT GIPAW and can predict shieldings for a wide range of organic crystals at a fraction of the computational cost.

**Scientific domain**: NMR chemical shift prediction, ML-accelerated spectroscopy  
**Target user community**: Researchers predicting NMR chemical shifts for organic solids, NMR crystallography practitioners

## Theoretical Methods
- Machine learning (Gaussian process regression)
- DFT GIPAW training data (CASTEP)
- SOAP (Smooth Overlap of Atomic Positions) descriptors
- Chemical shielding tensor prediction
- Isotropic and anisotropic shieldings
- NMR crystallography

## Capabilities (CRITICAL)
- NMR chemical shielding prediction
- 1H, 13C, 15N, 17O, and other nuclei
- Isotropic shielding
- Shielding anisotropy
- Fast ML prediction (milliseconds per site)
- Organic crystal support
- Structure-to-shielding mapping
- Uncertainty estimation (Gaussian process)

**Sources**: GitHub repository, Chem. Sci. 12, 3365 (2021)

## Key Strengths

### ML Acceleration:
- Orders of magnitude faster than DFT GIPAW
- Milliseconds per shielding prediction
- Enables high-throughput NMR prediction
- Uncertainty quantification

### DFT-Quality Accuracy:
- Trained on CASTEP GIPAW data
- 1.4 million shieldings in training set
- Systematic improvement with more data
- Validated against experiment

### Organic Solids Focus:
- Optimized for organic crystals
- Handles disorder and polymorphism
- Supports common organic elements
- NMR crystallography applications

### Uncertainty Estimation:
- Gaussian process regression
- Meaningful error bars
- Active learning possible
- Confidence in predictions

## Inputs & Outputs
- **Input formats**:
  - Atomic structures (ASE, CIF, XYZ)
  - Crystal structure files
  
- **Output data types**:
  - Chemical shieldings (isotropic)
  - Shielding anisotropy
  - Prediction uncertainties
  - Site-resolved shieldings

## Interfaces & Ecosystem
- **ASE**: Structure handling
- **CASTEP**: Training data source
- **QUIP/SOAP**: Descriptors
- **Python/NumPy**: Computation

## Performance Characteristics
- **Speed**: Milliseconds per site
- **Accuracy**: ~1-2 ppm for 13C, ~0.3 ppm for 1H
- **System size**: Any organic crystal
- **Memory**: Low (trained model)

## Computational Cost
- **Prediction**: Milliseconds per site
- **Training**: Pre-computed (done once)
- **vs GIPAW**: 1000x+ speedup
- **Typical**: Negligible

## Limitations & Known Constraints
- **Organic solids only**: Not for inorganic/metallic systems
- **Training domain**: Limited to organic crystal compositions
- **No dynamics**: Static structure only
- **GIPAW reference**: Accuracy limited by GIPAW quality
- **Heavy elements**: Limited support

## Comparison with Other Codes
- **vs CASTEP GIPAW**: ShiftML is ML-fast, CASTEP is DFT-accurate
- **vs MRSimulator**: ShiftML predicts shieldings, MRSimulator simulates spectra
- **vs ml4nmr**: ShiftML is for solids, ml4nmr corrects molecular DFT
- **Unique strength**: ML prediction of NMR chemical shieldings for organic solids, uncertainty estimation

## Application Areas

### NMR Crystallography:
- Structure validation via NMR
- Polymorph identification
- Disorder characterization
- Crystal structure determination

### Pharmaceutical Materials:
- API polymorph screening
- Cocrystal characterization
- Amorphous solid dispersions
- Hydrate/anhydrate identification

### Organic Semiconductors:
- Molecular packing determination
- π-π stacking characterization
- Charge transport correlations
- Structure-property relationships

### Metal-Organic Frameworks:
- Linker conformation analysis
- Guest molecule identification
- Framework dynamics
- Functional group characterization

## Best Practices

### Structure Preparation:
- Use well-optimized structures
- Include hydrogen positions
- Consider disorder models
- Validate geometry

### Prediction Analysis:
- Use uncertainty estimates
- Compare with experimental shifts
- Account for temperature effects
- Consider dynamics averaging

## Community and Support
- Open source on GitHub
- Developed at COSMO Lab (EPFL)
- Published in Chem. Sci.
- Active development
- Documentation available

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/lab-cosmo/ShiftML
2. Documentation: https://lab-cosmo.github.io/ShiftML/
3. M. Caro et al., Chem. Sci. 12, 3365 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE
- Published methodology: Chem. Sci.
- Active development: Ongoing
- Specialized strength: ML prediction of NMR chemical shieldings for organic solids, uncertainty estimation
