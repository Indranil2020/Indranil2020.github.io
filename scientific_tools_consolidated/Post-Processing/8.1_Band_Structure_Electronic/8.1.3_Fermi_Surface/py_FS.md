# py_FS

## Overview
**py_FS** is a Python script for plotting Fermi surfaces from BXSF files. It provides 3D visualization using PyVista.

## Official Resources
- **GitHub**: https://github.com/TheoWeinberger/py_FS

## Capabilities
- **Fermi Surface Plotting**: 3D visualization from BXSF
- **Multi-band Support**: Plot multiple bands
- **Interactive Visualization**: PyVista-based rendering

## Key Features
- BXSF file support
- PyVista 3D rendering
- Matplotlib integration
- Simple command-line usage

## Requirements
- Python 3.9+
- numpy, matplotlib, scipy
- pyvista, trimesh

## Usage
```bash
python fs_plot.py case
```
Where case.bxsf.band-n files are in the directory.

## Inputs & Outputs
- **Inputs**: BXSF files (XCrySDen/Quantum ESPRESSO format)
- **Outputs**: 3D Fermi surface visualizations

## Verification & Sources
- **Status**: ✅ VERIFIED
- **Confidence**: VERIFIED
