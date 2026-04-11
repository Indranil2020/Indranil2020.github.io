# CASCADE

## Official Resources
- Source Repository: https://github.com/patonlab/CASCADE
- Documentation: Included in repository
- License: Open source

## Overview
**CASCADE** (CAlculation of NMR Chemical Shifts using DEep learning) is a tool for predicting NMR chemical shifts using machine learning. It combines DFT-calculated shifts with ML corrections to achieve CCSD(T)-quality predictions at DFT cost, particularly for organic molecules.

**Scientific domain**: NMR chemical shift prediction, ML-accelerated spectroscopy  
**Target user community**: Researchers predicting NMR chemical shifts for organic molecules with ML-enhanced accuracy

## Theoretical Methods
- DFT NMR chemical shift calculation
- Machine learning correction (Δ-ML)
- CCSD(T) quality from DFT+ML
- 1H and 13C chemical shift prediction
- Spin-orbit relativistic corrections (ΔSO)

## Capabilities (CRITICAL)
- NMR chemical shift prediction
- ML correction of DFT shifts
- 1H and 13C chemical shifts
- CCSD(T)-quality accuracy at DFT cost
- Spin-orbit relativistic corrections
- Organic molecule support

**Sources**: GitHub repository, J. Chem. Inf. Model.

## Key Strengths

### ML-Enhanced Accuracy:
- CCSD(T)-quality predictions
- DFT cost with ML correction
- Systematic improvement
- Trained on high-level data

### NMR-Specific:
- 1H and 13C chemical shifts
- Spin-orbit corrections
- Heavy atom effects
- Organic molecule focus

### DFT Integration:
- Uses DFT as baseline
- ML correction on top
- Any DFT functional as input
- Flexible framework

## Inputs & Outputs
- **Input formats**:
  - Molecular geometry
  - DFT-calculated shifts
  - ML model files
  
- **Output data types**:
  - Corrected chemical shifts
  - ML correction values
  - Comparison with DFT
  - Prediction accuracy

## Interfaces & Ecosystem
- **DFT codes**: Baseline shift calculation
- **PyTorch/scikit-learn**: ML backend
- **Python**: Core language

## Performance Characteristics
- **Speed**: Fast (ML prediction)
- **Accuracy**: CCSD(T)-quality
- **System size**: Organic molecules
- **Memory**: Low

## Computational Cost
- **ML prediction**: Seconds
- **DFT baseline**: Minutes (separate)
- **Typical**: Very efficient

## Limitations & Known Constraints
- **Organic molecules**: Limited to organic systems
- **1H and 13C only**: No other nuclei
- **Training data dependent**: Quality limited by training set
- **DFT baseline needed**: Requires DFT calculation first

## Comparison with Other Codes
- **vs ShiftML**: CASCADE is organic-focused, ShiftML is solid-state NMR
- **vs ml4nmr**: CASCADE is Δ-ML correction, ml4nmr is direct ML correction
- **vs afnmr**: CASCADE is ML, afnmr is bio-NMR specific
- **Unique strength**: ML-corrected NMR chemical shifts to CCSD(T) quality for organic molecules

## Application Areas

### Organic Chemistry:
- Structure verification
- NMR shift prediction
- Isomer identification
- Reaction monitoring

### Drug Discovery:
- Small molecule NMR
- Metabolite identification
- Purity assessment
- Structural elucidation

### Method Development:
- ML-NMR benchmarking
- DFT functional comparison
- Chemical shift databases
- Accuracy improvement

## Best Practices

### DFT Baseline:
- Use consistent DFT functional
- Adequate basis set for NMR
- Include solvent effects
- Validate against experiment

### ML Correction:
- Use appropriate ML model
- Check domain of applicability
- Validate on test set
- Compare DFT vs corrected

## Community and Support
- Open source on GitHub
- Developed by Paton Lab (CSU)
- Published in J. Chem. Inf. Model.
- Research code

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/patonlab/CASCADE

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Published methodology: J. Chem. Inf. Model.
- Specialized strength: ML-corrected NMR chemical shifts to CCSD(T) quality for organic molecules
