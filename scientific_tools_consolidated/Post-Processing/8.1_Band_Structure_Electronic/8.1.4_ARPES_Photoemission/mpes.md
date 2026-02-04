# mpes

## Official Resources
- **Homepage**: https://github.com/mpes-kit/mpes
- **GitHub**: https://github.com/mpes-kit/mpes
- **Documentation**: https://mpes-kit.github.io/mpes/
- **PyPI**: https://pypi.org/project/mpes/
- **License**: MIT License

## Overview
mpes is a Python toolkit for multidimensional photoemission spectroscopy data analysis. It provides comprehensive tools for processing, calibrating, and analyzing time-resolved and angle-resolved photoemission data, with support for modern data formats and parallel processing.

**Scientific domain**: ARPES data analysis, time-resolved photoemission, ultrafast dynamics
**Target user community**: Experimental physicists working with ARPES and trARPES data

## Theoretical Background
mpes processes photoemission data based on:
- Angle-to-momentum conversion: k = (1/ℏ)√(2mE_kin) sin(θ)
- Energy calibration using Fermi edge fitting
- Time-resolved dynamics analysis
- Multidimensional data binning and interpolation

## Capabilities (CRITICAL)
- **Data Processing**: ARPES/trARPES data handling and binning
- **Calibration**: Energy, momentum, and delay calibration
- **Visualization**: Multi-dimensional data plotting
- **Analysis**: Band structure and dynamics extraction
- **File I/O**: HDF5, FITS, and beamline-specific formats
- **Parallel Processing**: Dask integration for large datasets

## Key Strengths

### Multidimensional Support:
- 2D, 3D, 4D data handling
- Time-resolved ARPES (trARPES)
- Pump-probe delay scans
- Temperature-dependent measurements

### Calibration Tools:
- Fermi edge fitting
- Momentum calibration
- Delay zero determination
- Distortion correction

### Modern Data Handling:
- HDF5 native support
- Dask for parallel processing
- xarray integration
- Memory-efficient binning

## Inputs & Outputs
- **Input formats**:
  - HDF5 files
  - FITS files
  - Beamline-specific formats
  - Raw detector data
  
- **Output data types**:
  - Calibrated spectra
  - Momentum-resolved data
  - Time-resolved dynamics
  - Matplotlib/Plotly figures

## Installation
```bash
pip install mpes
```

## Usage Examples
```python
import mpes

# Load data
data = mpes.io.load_hdf5("arpes_data.h5")

# Calibrate energy
calibrated = mpes.calibration.calibrate_energy(data, fermi_level=16.8)

# Convert to momentum
kdata = mpes.conversion.angle_to_k(calibrated)

# Plot
mpes.visualization.plot_spectrum(kdata)
```

## Performance Characteristics
- **Speed**: Dask parallel processing
- **Memory**: Efficient for large datasets
- **Scalability**: Handles multi-GB datasets

## Limitations & Known Constraints
- **Experimental focus**: Designed for experimental data
- **Beamline-specific**: Some features beamline-dependent
- **Learning curve**: Complex calibration procedures

## Comparison with Other Tools
- **vs PyARPES**: Different feature sets, both comprehensive
- **vs ARPESGUI**: mpes is Python-native, scriptable
- **Unique strength**: trARPES support, Dask integration

## Application Areas
- ARPES experiments
- Time-resolved photoemission (trARPES)
- Ultrafast dynamics studies
- Band structure mapping
- Fermi surface measurements

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mpes-kit/mpes

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Developer: mpes-kit team
