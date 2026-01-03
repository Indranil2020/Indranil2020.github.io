# PyARPES

## Official Resources
- Homepage: https://github.com/chstan/arpes
- Documentation: https://arpes.readthedocs.io/
- Source Repository: https://github.com/chstan/arpes
- License: MIT License

## Overview
PyARPES is a Python-based analysis framework for Angle-Resolved Photoemission Spectroscopy (ARPES) data. It provides tools for loading, processing, visualizing, and analyzing multidimensional ARPES datasets. PyARPES aims to streamline the workflow from raw data to publication-quality figures, supporting data from various synchrotrons and lab-based systems.

**Scientific domain**: ARPES data analysis, photoemission spectroscopy, electronic structure  
**Target user community**: Experimental condensed matter physicists, ARPES practitioners

## Capabilities (CRITICAL)
- **Data Loading**: Support for various file formats (HDF5, Igor, Fits, etc.) from major beamlines (ALS, SSRL, Diamond, etc.)
- **Multidimensional Analysis**: Handling of 2D, 3D, and 4D datasets (Energy, kx, ky, time/temperature)
- **Visualization**: Interactive and publication-quality plotting
- **Processing**: Normalization, background subtraction, curvature analysis
- **Curve Fitting**: MDC/EDC fitting, Fermi edge fitting
- **Coordinate Transformation**: Conversion between angles and momentum space
- **Automation**: Scriptable workflows for batch processing

**Sources**: PyARPES documentation, GitHub repository

## Inputs & Outputs
- **Input formats**: .ibw (Igor Binary), .h5/.nxs (HDF5/NeXus), .fits, .txt, .zip
- **Output data types**: Processed datasets (HDF5), Matplotlib figures, fitted parameters

## Interfaces & Ecosystem
- **Python**: Built on xarray, pandas, and matplotlib
- **Jupyter**: Designed for use in Jupyter notebooks
- **xarray**: Uses xarray for labeled multidimensional arrays

## Workflow and Usage
1. Load data: `data = load_data("file.ibw")`
2. Preprocess: `data = data.arpes.k_convert()` (convert angles to k-space)
3. Visualize: `data.sum("kx").plot()`
4. Fit: Extract dispersions or gaps using built-in fitting tools

## Performance Characteristics
- Python-based, leverages numpy/pandas for efficiency
- Handles large datasets via lazy loading (dask integration possible)

## Application Areas
- Band structure mapping
- Fermi surface mapping
- Gap analysis in superconductors
- Time-resolved ARPES analysis

## Community and Support
- Open-source (MIT)
- GitHub issues for support
- Developed by researchers at Stanford/UBC

## Verification & Sources
**Primary sources**:
1. Homepage: https://github.com/chstan/arpes
2. Documentation: https://arpes.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Research community)
- Applications: ARPES data analysis, visualization, coordinate transformation
