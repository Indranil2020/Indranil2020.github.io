# ml4nmr

## Official Resources
- Source Repository: https://github.com/grimme-lab/ml4nmr
- Documentation: Included in repository
- License: Open source

## Overview
**ml4nmr** is a machine learning-based correction tool for NMR chemical shifts calculated with DFT. It enables correction of 1H and 13C NMR chemical shifts toward CCSD(T) quality (ΔcorrML) and prediction of spin-orbit relativistic contributions to NMR chemical shifts caused by heavy atoms (ΔSO^ML).

**Scientific domain**: ML-corrected NMR chemical shifts, spin-orbit corrections  
**Target user community**: Researchers needing high-accuracy NMR chemical shift predictions including heavy-atom relativistic effects

## Theoretical Methods
- Machine learning correction of DFT NMR shifts (ΔcorrML)
- CCSD(T) quality from DFT+ML
- Spin-orbit relativistic correction (ΔSO^ML)
- 1H and 13C chemical shift correction
- Heavy atom effects on chemical shifts
- Grimme lab DFT baseline (XTB)

## Capabilities (CRITICAL)
- ML correction of DFT NMR chemical shifts
- ΔcorrML: DFT→CCSD(T) quality correction
- ΔSO^ML: Spin-orbit relativistic correction
- 1H and 13C chemical shifts
- Heavy atom relativistic effects
- XTB/DFT integration

**Sources**: GitHub repository, J. Chem. Theory Comput.

## Key Strengths

### Dual ML Correction:
- ΔcorrML: Accuracy correction to CCSD(T)
- ΔSO^ML: Relativistic correction for heavy atoms
- Combined: High accuracy + relativistic effects
- Systematic improvement

### Heavy Atom Support:
- Spin-orbit corrections
- Relativistic effects on shifts
- Halogen and heavy element effects
- Beyond DFT accuracy

### Grimme Lab Quality:
- Well-tested ML models
- XTB integration
- Consistent with Grimme ecosystem
- Published benchmarks

## Inputs & Outputs
- **Input formats**:
  - DFT/XTB calculated shifts
  - Molecular geometry
  - ML model files
  
- **Output data types**:
  - Corrected chemical shifts
  - ΔcorrML correction values
  - ΔSO^ML correction values
  - Final corrected shifts

## Interfaces & Ecosystem
- **XTB**: Grimme's semi-empirical code
- **DFT codes**: Baseline shift calculation
- **Python/PyTorch**: ML backend

## Performance Characteristics
- **Speed**: Fast (ML prediction)
- **Accuracy**: CCSD(T) + relativistic quality
- **System size**: Organic molecules
- **Memory**: Low

## Computational Cost
- **ML prediction**: Seconds
- **DFT/XTB baseline**: Minutes (separate)
- **Typical**: Very efficient

## Limitations & Known Constraints
- **1H and 13C only**: No other nuclei
- **Organic focus**: Limited inorganic support
- **Training dependent**: Quality limited by training data
- **XTB preferred**: Best with Grimme ecosystem

## Comparison with Other Codes
- **vs CASCADE**: ml4nmr has ΔSO^ML relativistic, CASCADE is Δ-ML accuracy
- **vs ShiftML**: ml4nmr is molecular, ShiftML is solid-state
- **vs afnmr**: ml4nmr is ML, afnmr is bio-NMR
- **Unique strength**: ML correction of NMR shifts with spin-orbit relativistic effects (ΔSO^ML)

## Application Areas

### Organic Chemistry:
- High-accuracy NMR prediction
- Heavy atom effects on shifts
- Structure verification
- Isomer identification

### Organometallic Chemistry:
- Heavy atom relativistic effects
- Metal-containing molecules
- Halogen effects on shifts
- Spin-orbit coupling

### Method Development:
- ML-NMR benchmarking
- Relativistic correction validation
- Chemical shift databases
- DFT functional assessment

## Best Practices

### DFT Baseline:
- Use consistent functional
- Include all atoms in calculation
- Use appropriate basis set
- Validate against experiment

### ML Correction:
- Apply both ΔcorrML and ΔSO^ML
- Check domain of applicability
- Validate on test set
- Compare corrected vs uncorrected

## Community and Support
- Open source on GitHub
- Developed by Grimme Lab (Bonn)
- Published in J. Chem. Theory Comput.
- Research code

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/grimme-lab/ml4nmr

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Published methodology: J. Chem. Theory Comput.
- Specialized strength: ML correction of NMR shifts with spin-orbit relativistic effects (ΔSO^ML)
