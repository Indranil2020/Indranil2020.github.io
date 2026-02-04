# IFermi

## Overview
**IFermi** is a Python package for Fermi surface generation, analysis, and visualization. It provides tools for extracting Fermi surfaces from DFT calculations and creating publication-quality visualizations.

## Official Resources
- **GitHub**: https://github.com/fermisurfaces/IFermi
- **Documentation**: https://fermisurfaces.github.io/IFermi/
- **PyPI**: https://pypi.org/project/ifermi/

## Capabilities
- **Fermi Surface Generation**: Extract from DFT band structures
- **2D Slices**: Fermi surface cross-sections
- **3D Visualization**: Interactive Fermi surface plots
- **Spin Textures**: Spin-resolved Fermi surfaces
- **Analysis**: Fermi surface properties calculation

## Key Features
- VASP, Quantum ESPRESSO support
- Interactive visualization (mayavi, plotly, matplotlib)
- Spin texture plotting
- Command-line interface
- Python API

## Supported Codes
- VASP (vasprun.xml)
- Quantum ESPRESSO (via BoltzTraP2)
- Any code via BXSF format

## Installation
```bash
pip install ifermi
```

## Usage
```bash
ifermi plot vasprun.xml
ifermi plot --slice 0 0 1 vasprun.xml
```

## Verification & Sources
- **Status**: ✅ VERIFIED
- **Confidence**: VERIFIED
- **Developer**: Materials Project team
