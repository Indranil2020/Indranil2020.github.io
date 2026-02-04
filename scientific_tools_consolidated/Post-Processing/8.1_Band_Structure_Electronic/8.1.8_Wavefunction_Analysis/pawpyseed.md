# pawpyseed

## Official Resources
- **Homepage**: https://github.com/kylebystrom/pawpyseed
- **GitHub**: https://github.com/kylebystrom/pawpyseed
- **Documentation**: https://pawpyseed.readthedocs.io/
- **PyPI**: https://pypi.org/project/pawpyseed/
- **License**: BSD 3-Clause License

## Overview
pawpyseed is a parallel C/Python package for numerical analysis of PAW (Projector Augmented Wave) DFT wavefunctions from VASP calculations. It enables reconstruction of all-electron wavefunctions from pseudo-wavefunctions and PAW projectors, essential for accurate defect analysis and wavefunction overlap calculations.

**Scientific domain**: Wavefunction analysis, PAW method, defect physics
**Target user community**: Researchers studying defects, core-level properties, and wavefunction character in VASP calculations

## Theoretical Background
pawpyseed implements PAW wavefunction reconstruction:
- All-electron wavefunction: |ψ_AE⟩ = |ψ̃⟩ + Σ_i (|φ_i⟩ - |φ̃_i⟩)⟨p̃_i|ψ̃⟩
- Pseudo-wavefunction from WAVECAR
- PAW projectors from POTCAR
- Core reconstruction for accurate overlaps
- Defect wavefunction localization analysis

## Capabilities (CRITICAL)
- **Wavefunction Reconstruction**: Full all-electron wavefunctions from PAW
- **PAW Analysis**: Projector function handling and core reconstruction
- **Overlap Calculations**: Accurate wavefunction overlaps including core
- **Defect Analysis**: Defect wavefunction characterization and localization
- **Band Decomposition**: Decompose defect states into bulk bands
- **Charge Analysis**: Core and valence charge separation
- **Parallel Processing**: C backend for performance

## Key Strengths

### PAW Reconstruction:
- Full all-electron wavefunctions
- Core state reconstruction
- Accurate for core-level properties
- Proper PAW augmentation

### Defect Physics:
- Defect wavefunction localization
- Bulk band decomposition
- Transition level analysis
- Configuration coordinate diagrams

### Performance:
- C backend for speed
- Parallel processing
- Efficient WAVECAR reading
- Memory-optimized

### Python Interface:
- Clean Python API
- NumPy integration
- Scriptable analysis
- Jupyter compatible

## Inputs & Outputs
- **Input formats**:
  - WAVECAR (pseudo-wavefunctions)
  - POTCAR (PAW projectors)
  - POSCAR (structure)
  
- **Output data types**:
  - All-electron wavefunctions
  - Overlap matrices
  - Charge densities
  - Decomposition coefficients

## Interfaces & Ecosystem
- **Python integration**:
  - NumPy for arrays
  - SciPy for linear algebra
  - Matplotlib for visualization
  
- **VASP compatibility**:
  - Standard VASP output files
  - PAW pseudopotentials
  - Gamma and k-point calculations

## Installation
```bash
pip install pawpyseed
```

From source (for development):
```bash
git clone https://github.com/kylebystrom/pawpyseed.git
cd pawpyseed
pip install -e .
```

## Usage Examples
```python
from pawpyseed.core.wavefunction import Wavefunction
from pawpyseed.core.projector import Projector

# Load wavefunction
wf = Wavefunction.from_directory("path/to/vasp/calc")

# Get all-electron wavefunction
ae_wf = wf.get_ae_wavefunction(band=0, kpoint=0, spin=0)

# Calculate overlap
overlap = wf.band_overlap(wf2, band1=0, band2=0)

# Defect analysis
from pawpyseed.analysis.defect import DefectAnalysis
defect = DefectAnalysis(bulk_wf, defect_wf)
decomposition = defect.get_decomposition()
```

## Performance Characteristics
- **Speed**: C backend, parallel processing
- **Memory**: Efficient for large WAVECAR files
- **Accuracy**: Full PAW reconstruction
- **Scalability**: Handles large supercells

## Limitations & Known Constraints
- **VASP-specific**: Only works with VASP PAW calculations
- **POTCAR required**: Needs PAW projector data
- **Memory**: Large supercells need significant RAM
- **Compilation**: C extension requires compiler

## Comparison with Other Tools
- **vs pymatgen**: pawpyseed specialized for PAW wavefunctions
- **vs VaspBandUnfolding**: Different focus (PAW vs unfolding)
- **Unique strength**: Full all-electron PAW reconstruction, defect analysis

## Application Areas
- Point defect characterization
- Defect transition levels
- Core-level spectroscopy
- Wavefunction localization
- Band decomposition analysis
- Configuration coordinate diagrams

## Best Practices
- Use consistent POTCAR for bulk and defect
- Check wavefunction normalization
- Verify overlap convergence
- Use appropriate supercell size

## Community and Support
- GitHub issue tracker
- Documentation with examples
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/kylebystrom/pawpyseed
2. K. Bystrom et al., related publications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (GitHub, BSD 3-Clause)
- Developer: Kyle Bystrom
- Active development: Maintained
