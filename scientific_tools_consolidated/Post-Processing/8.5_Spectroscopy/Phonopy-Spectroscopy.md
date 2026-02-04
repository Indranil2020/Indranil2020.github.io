# Phonopy-Spectroscopy

## Official Resources
- Homepage: https://github.com/skelton-group/Phonopy-Spectroscopy
- GitHub: https://github.com/skelton-group/Phonopy-Spectroscopy
- Documentation: GitHub README and examples
- License: MIT License

## Overview
Phonopy-Spectroscopy is a Python package that extends Phonopy to simulate infrared (IR) and Raman spectra of crystalline materials. It calculates IR intensities from Born effective charges and Raman activities from polarizability tensors computed using VASP, enabling direct comparison with experimental vibrational spectra.

**Scientific domain**: IR spectroscopy, Raman spectroscopy, phonons
**Target user community**: Materials scientists studying vibrational properties

## Theoretical Methods
- Phonon frequencies (via Phonopy)
- Born effective charge tensors
- Dielectric polarizability tensors
- IR intensity calculation
- Raman activity calculation
- Symmetry analysis

## Capabilities (CRITICAL)
- **IR Spectra**: Infrared intensities from Born charges
- **Raman Spectra**: Raman activities from polarizabilities
- **Phonopy Integration**: Uses Phonopy phonon modes
- **VASP Interface**: Reads VASP output directly
- **Peak Broadening**: Lorentzian/Gaussian functions
- **Symmetry**: Irreducible representation assignment

**Sources**: Phonopy-Spectroscopy GitHub, Skelton group publications

## Key Strengths

### Phonopy Ecosystem:
- Built on Phonopy
- Consistent workflow
- Validated methodology
- Active maintenance

### VASP Workflow:
- Direct VASP integration
- Born charges (LEPSILON)
- Raman from DFPT
- Standard DFT inputs

### Complete Spectra:
- IR and Raman
- Peak assignment
- Symmetry labels
- Publication ready

## Inputs & Outputs
- **Input formats**:
  - Phonopy FORCE_SETS
  - VASP OUTCAR (Born charges)
  - VASP polarizability calculations
  
- **Output data types**:
  - IR spectrum data
  - Raman spectrum data
  - Mode assignments
  - Visualization files

## Installation
```bash
pip install phonopy-spectroscopy
# Or from source
git clone https://github.com/skelton-group/Phonopy-Spectroscopy.git
cd Phonopy-Spectroscopy
pip install -e .
```

## Usage Examples
```bash
# Calculate IR spectrum
phonopy-ir POSCAR --born BORN --fc FORCE_CONSTANTS

# Calculate Raman spectrum
phonopy-raman POSCAR --raman-tensors tensors.yaml --fc FORCE_CONSTANTS

# Python API
from spectroscopy import calculate_ir_spectrum
spectrum = calculate_ir_spectrum(phonopy_obj, born_charges)
```

## Performance Characteristics
- **Speed**: Fast post-processing
- **Memory**: Minimal requirements
- **Dependencies**: Requires Phonopy and DFT data

## Limitations & Known Constraints
- **VASP focused**: Primarily for VASP users
- **Phonopy required**: Needs Phonopy installation
- **DFPT needed**: Requires response calculations
- **Harmonic**: No anharmonic effects

## Comparison with Other Tools
- **vs THeSeuSS**: Different implementations
- **vs native DFT**: More automated workflow
- **vs experiment**: Direct comparison enabled
- **Unique strength**: Phonopy integration, VASP workflow

## Application Areas
- Materials characterization
- Phase identification
- Catalyst analysis
- Mineral spectroscopy
- Polymer identification

## Best Practices
- Ensure phonon convergence first
- Use LEPSILON for Born charges
- Check symmetry assignments
- Compare with experimental spectra

## Community and Support
- GitHub repository
- Skelton group (University of Manchester)
- MIT licensed
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/skelton-group/Phonopy-Spectroscopy
2. Skelton group publications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (MIT)
- Method: IR/Raman from Phonopy
