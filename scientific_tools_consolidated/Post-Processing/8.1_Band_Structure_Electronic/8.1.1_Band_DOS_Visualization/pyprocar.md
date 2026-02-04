# pyprocar

## Official Resources
- **Homepage**: https://romerogroup.github.io/pyprocar/
- **GitHub**: https://github.com/romerogroup/pyprocar
- **Documentation**: https://pyprocar.readthedocs.io/
- **PyPI**: https://pypi.org/project/pyprocar/
- **License**: GPL v3

## Overview
pyprocar is a comprehensive Python library for pre- and post-processing of electronic structure data from multiple DFT codes. It excels at plotting band structures with orbital projections, Fermi surfaces, and spin textures with publication-quality output.

**Scientific domain**: Electronic structure visualization, multi-code post-processing
**Target user community**: DFT users (VASP, QE, Elk, ABINIT) needing advanced visualization

## Theoretical Background
pyprocar processes electronic structure data for:
- Projected band structures (fatbands)
- Fermi surface visualization
- Spin texture analysis
- Band unfolding from supercells

## Capabilities (CRITICAL)
- **Band Structure**: Standard and projected band plotting
- **Fermi Surface**: 2D/3D Fermi surface visualization
- **Spin Textures**: Spin-resolved band analysis
- **Band Unfolding**: Supercell to primitive cell unfolding
- **DOS/PDOS**: Density of states analysis
- **Multi-Code**: VASP, QE, Elk, ABINIT support
- **Publication Quality**: High-quality matplotlib figures

## Key Strengths

### Multi-Code Support:
- VASP (PROCAR, vasprun.xml)
- Quantum ESPRESSO
- Elk
- ABINIT

### Advanced Visualization:
- Orbital-projected fatbands
- Spin texture mapping
- 3D Fermi surfaces
- Customizable color schemes

### Band Unfolding:
- Supercell unfolding
- Spectral weight visualization
- Defect/alloy band structures

## Inputs & Outputs
- **Input formats**: PROCAR, vasprun.xml (VASP), QE output, Elk output, ABINIT output
- **Output data types**: Matplotlib figures, publication-ready plots

## Installation
```bash
pip install pyprocar
```

## Usage Examples
```python
import pyprocar

# Plot band structure
pyprocar.bandsplot('PROCAR', outcar='OUTCAR', mode='plain')

# Projected bands (fatbands)
pyprocar.bandsplot('PROCAR', outcar='OUTCAR', mode='parametric',
                   orbitals=[0,1,2], atoms=[0])

# Fermi surface
pyprocar.fermi3D('PROCAR', outcar='OUTCAR', bands=[4,5])

# Band unfolding
pyprocar.unfold('PROCAR', outcar='OUTCAR', supercell_matrix=[[2,0,0],[0,2,0],[0,0,2]])
```

## Performance Characteristics
- **Speed**: Efficient Python implementation
- **Memory**: Handles large PROCAR files
- **Visualization**: High-quality matplotlib output

## Limitations & Known Constraints
- **PROCAR required**: Needs projection data
- **Memory**: Large supercells need significant RAM
- **Learning curve**: Advanced features require practice

## Comparison with Other Tools
- **vs sumo**: pyprocar has more projection features
- **vs vaspvis**: pyprocar multi-code, vaspvis VASP-only
- **Unique strength**: Multi-code support, comprehensive projections

## Application Areas
- Orbital character analysis
- Spin texture visualization
- Fermi surface studies
- Defect band structures
- Alloy electronic structure

## Verification & Sources
**Primary sources**:
1. Documentation: https://pyprocar.readthedocs.io/
2. GitHub: https://github.com/romerogroup/pyprocar

**Confidence**: VERIFIED - Developed by Romero Group

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, GPL v3)
- Developer: Romero Group (WVU)
- Active development: Regular releases
