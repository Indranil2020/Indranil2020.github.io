# Larch (xraylarch)

## Official Resources
- Homepage: https://xraypy.github.io/xraylarch/
- GitHub: https://github.com/xraypy/xraylarch
- Documentation: https://xraypy.github.io/xraylarch/
- PyPI: https://pypi.org/project/xraylarch/
- Publication: M. Newville, J. Synchrotron Rad. 20, 913 (2013)
- License: BSD 3-Clause License

## Overview
Larch is an open-source Python library and set of applications for processing and analyzing X-ray absorption spectroscopy (XAS) data from synchrotron beamlines. It provides comprehensive tools for XAFS (X-ray Absorption Fine Structure), including XANES (near-edge) and EXAFS (extended) analysis, as well as XRF (X-ray fluorescence) mapping and analysis.

**Scientific domain**: X-ray absorption spectroscopy, XAFS, XANES, EXAFS, XRF
**Target user community**: Synchrotron beamline scientists and XAS researchers

## Theoretical Methods
- XAFS data processing and analysis
- Background subtraction and normalization
- Fourier transform for EXAFS
- FEFF path fitting
- Linear combination fitting
- Principal component analysis
- XRF peak fitting

## Capabilities (CRITICAL)
- **XANES Analysis**: Pre-edge, edge, and post-edge analysis
- **EXAFS Analysis**: Fourier transform, fitting with FEFF paths
- **XRF Mapping**: Fluorescence data processing
- **GUI Applications**: Larix (XAS Viewer), GSE MapViewer
- **FEFF Integration**: Path fitting with FEFF calculations
- **Batch Processing**: Scripting for high-throughput
- **Database**: X-ray absorption edge and emission line data

**Sources**: Larch documentation, J. Synchrotron Rad. publication

## Key Strengths

### Comprehensive XAS:
- Full XAFS workflow
- XANES and EXAFS analysis
- Background removal
- Normalization and merging

### Python Native:
- Modern Python implementation
- NumPy/SciPy based
- Jupyter compatible
- Scriptable workflows

### GUI Applications:
- Larix for XAS viewing
- Interactive fitting
- Visualization tools
- User-friendly interface

## Inputs & Outputs
- **Input formats**:
  - ASCII data files
  - HDF5 files
  - Various beamline formats
  - FEFF output files
  
- **Output data types**:
  - Processed spectra
  - Fit results
  - χ(k) and χ(R) data
  - Reports and figures

## Installation
```bash
pip install xraylarch
# Or with conda
conda install -c conda-forge xraylarch
```

## Usage Examples
```python
from larch import Interpreter
from larch.xafs import autobk, xftf

# Create Larch session
session = Interpreter()

# Read and process XAS data
data = read_ascii('fe_foil.dat')

# Background subtraction
autobk(data, rbkg=1.0, kweight=2)

# Fourier transform
xftf(data, kmin=2, kmax=12, dk=2, window='hanning')

# Plot results
plot_chir(data)
```

## Performance Characteristics
- **Speed**: Efficient Python/NumPy implementation
- **Memory**: Handles large datasets
- **Parallelization**: Batch processing support

## Limitations & Known Constraints
- **XAS focus**: Primarily for X-ray absorption
- **Learning curve**: XAS concepts required
- **FEFF needed**: Path fitting requires FEFF calculations
- **Beamline formats**: Some formats need custom readers

## Comparison with Other Tools
- **vs Demeter**: Larch Python-native, Demeter Perl-based
- **vs FEFF**: Larch analysis, FEFF calculation
- **vs Athena/Artemis**: Larch is modern successor
- **Unique strength**: Python ecosystem, comprehensive XAS

## Application Areas
- Synchrotron XAS experiments
- Catalysis research
- Environmental science
- Materials characterization
- Coordination chemistry

## Best Practices
- Calibrate energy scale
- Use appropriate background parameters
- Validate with reference compounds
- Document processing parameters

## Community and Support
- GitHub repository
- Mailing list (XAFS community)
- Active development
- Matt Newville (developer)

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/xraypy/xraylarch
2. M. Newville, J. Synchrotron Rad. 20, 913 (2013)
3. Documentation: https://xraypy.github.io/xraylarch/

**Confidence**: VERIFIED - Published in J. Synchrotron Rad.

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (BSD-3)
- Academic citations: Well-cited
- Active development: Maintained
