# XANESNET

## Official Resources
- Source Repository: https://github.com/NewcastleRSE/xray-spectroscopy-ml
- License: Open source

## Overview
**XANESNET** is a deep neural network for predicting X-ray Absorption Near Edge Structure (XANES) spectra from molecular structure descriptors. It enables fast and accurate prediction of transition metal XAS spectra without requiring expensive first-principles calculations, making it suitable for high-throughput screening and real-time spectral prediction.

**Scientific domain**: ML-accelerated X-ray spectroscopy, XANES prediction  
**Target user community**: Researchers needing rapid XANES prediction for transition metal systems, ML-spectroscopy developers

## Theoretical Methods
- Deep neural network (DNN) for spectra prediction
- Molecular structure descriptors (FEEF, SOAP, etc.)
- XANES spectral prediction
- Valence-to-core XES prediction
- Transfer learning approaches
- Regression on spectral features

## Capabilities (CRITICAL)
- K-edge XANES prediction for transition metals
- L-edge XANES prediction
- Valence-to-core XES prediction
- Structure-to-spectrum mapping
- High-throughput spectral screening
- Real-time spectral prediction
- Training on FEFF-calculated or experimental data

**Sources**: J. Chem. Phys. 156, 164102 (2022)

## Key Strengths

### Speed:
- Milliseconds per spectrum prediction
- Orders of magnitude faster than FEFF/DFT
- Enables real-time analysis
- High-throughput screening feasible

### Accuracy:
- Trained on FEFF or experimental data
- Good generalization within training domain
- Systematic improvement with more data
- Competitive with first-principles for standard edges

### Flexibility:
- Multiple edge types
- Multiple descriptor types
- Transfer learning between edges
- Custom training data

## Inputs & Outputs
- **Input formats**:
  - Molecular structure descriptors
  - XYZ coordinates
  - Pre-computed FEFF descriptors
  
- **Output data types**:
  - Predicted XANES spectra
  - Predicted VtC XES spectra
  - Confidence metrics
  - Feature importance analysis

## Interfaces & Ecosystem
- **FEFF**: Training data generation
- **Python/PyTorch**: ML framework
- **ASE**: Structure handling
- **NumPy**: Data processing

## Performance Characteristics
- **Speed**: Milliseconds per prediction
- **Accuracy**: Good within training domain
- **System size**: Limited by descriptor, not prediction
- **Memory**: Low (trained model)

## Computational Cost
- **Prediction**: Milliseconds
- **Training**: Hours to days on GPU
- **FEFF training data**: Minutes per spectrum
- **Typical**: Very efficient after training

## Limitations & Known Constraints
- **Training data dependent**: Quality limited by training set
- **Extrapolation**: Poor outside training domain
- **No physics guarantee**: ML is interpolative
- **Transition metals**: Primarily validated for TM K-edges
- **Interpretability**: Limited (black box)

## Comparison with Other Codes
- **vs FEFF**: XANESNET is much faster, FEFF is more general
- **vs MLXANES**: XANESNET uses DNN, MLXANES uses linear regression
- **vs pyFitIt**: XANESNET is pure prediction, pyFitIt is fitting
- **Unique strength**: Deep neural network XANES prediction, fast and accurate for transition metals

## Application Areas

### High-Throughput Screening:
- Materials databases
- Catalyst screening
- Battery material discovery
- Rapid spectral assessment

### Real-Time Analysis:
- In situ/operando spectroscopy
- Beamline data interpretation
- Quick structure validation
- On-the-fly spectral prediction

### Inverse Problems:
- Spectrum-to-structure mapping
- Structural determination from XANES
- Fitting experimental spectra
- Constraint generation for refinement

## Best Practices

### Training Data:
- Use diverse structural motifs
- Include relevant oxidation states
- Validate on held-out test set
- Monitor overfitting

### Prediction:
- Check if input is within training domain
- Compare with FEFF for validation
- Use ensemble predictions for uncertainty
- Report confidence metrics

## Community and Support
- Open source on GitHub
- Developed at Newcastle University
- Published in J. Chem. Phys.
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/NewcastleRSE/xray-spectroscopy-ml
2. C. D. Rankine and T. J. Penfold, J. Chem. Phys. 156, 164102 (2022)
3. T. I. Madan et al., J. Chem. Phys. 159, 164102 (2023)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Published methodology: J. Chem. Phys.
- Active development: Ongoing
- Specialized strength: Deep neural network for XANES prediction, fast ML-accelerated spectroscopy
