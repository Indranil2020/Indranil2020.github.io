# HyperSpy

## Official Resources
- Homepage: https://hyperspy.org/
- GitHub: https://github.com/hyperspy/hyperspy
- Documentation: https://hyperspy.org/hyperspy-doc/current/
- Publication: F. de la Peña et al., Microsc. Microanal. 23, 214 (2017)
- License: GNU General Public License v3.0

## Overview
HyperSpy is an open-source Python library for multi-dimensional data analysis, with particular focus on electron microscopy data including EELS (Electron Energy Loss Spectroscopy), EDS (Energy Dispersive X-ray Spectroscopy), and other spectroscopic imaging techniques. It provides tools for data visualization, processing, and quantitative analysis.

**Scientific domain**: EELS, EDS, spectroscopic imaging analysis
**Target user community**: Electron microscopists and spectroscopists

## Theoretical Methods
- EELS quantification
- EDS quantification
- Principal component analysis (PCA)
- Independent component analysis (ICA)
- Machine learning decomposition
- Model fitting and curve fitting

## Capabilities (CRITICAL)
- **EELS Analysis**: Core-loss and low-loss processing
- **EDS Analysis**: Elemental quantification
- **Spectrum Imaging**: Multi-dimensional datasets
- **Decomposition**: PCA, ICA, NMF
- **Model Fitting**: Peak fitting, background removal
- **Visualization**: Interactive plotting
- **Big Data**: Lazy loading for large datasets

**Sources**: HyperSpy documentation, Microsc. Microanal. publication

## Key Strengths

### Multi-Dimensional:
- Spectrum images
- 4D-STEM data
- Time series
- Arbitrary dimensions

### Analysis Tools:
- Decomposition methods
- Curve fitting
- Elemental mapping
- Quantification

### Python Ecosystem:
- NumPy/SciPy based
- Jupyter compatible
- scikit-learn integration
- Active development

## Inputs & Outputs
- **Input formats**:
  - Digital Micrograph (.dm3, .dm4)
  - EMSA/MSA
  - HDF5, NetCDF
  - Many microscopy formats
  
- **Output data types**:
  - HyperSpy signals
  - HDF5 files
  - Matplotlib figures
  - Quantification results

## Installation
```bash
pip install hyperspy
# Or with conda
conda install -c conda-forge hyperspy
```

## Usage Examples
```python
import hyperspy.api as hs

# Load EELS spectrum image
s = hs.load('eels_data.dm4')

# Remove background
s.remove_background(signal_range=(500., 550.))

# Elemental quantification
s.quantification(method='CL')

# PCA decomposition
s.decomposition()
s.plot_decomposition_results()

# Fit a model
m = s.create_model()
m.fit()
```

## Performance Characteristics
- **Speed**: Efficient NumPy operations
- **Memory**: Lazy loading for big data
- **Scalability**: Handles large datasets

## Limitations & Known Constraints
- **Learning curve**: Many features to learn
- **TEM focus**: Primarily for electron microscopy
- **Dependencies**: Many package dependencies
- **Documentation**: Extensive but complex

## Comparison with Other Tools
- **vs Digital Micrograph**: HyperSpy open-source, scriptable
- **vs EELSMODEL**: Different approaches
- **vs NeXus/silx**: Different focus
- **Unique strength**: Comprehensive Python EELS/EDS

## Application Areas
- Electron microscopy
- Materials characterization
- Catalyst analysis
- Semiconductor analysis
- Thin film analysis

## Best Practices
- Use appropriate background models
- Validate quantification
- Consider artifacts
- Document processing steps

## Community and Support
- GitHub repository
- Gitter chat
- Mailing list
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/hyperspy/hyperspy
2. F. de la Peña et al., Microsc. Microanal. 23, 214 (2017)
3. Documentation: https://hyperspy.org/

**Confidence**: VERIFIED - Standard EELS analysis tool

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GPL-3.0)
- Academic citations: Well-cited
- Active development: Maintained
