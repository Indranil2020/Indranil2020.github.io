# vaspvis

## Official Resources
- **Homepage**: https://derekdardzinski.github.io/vaspvis/
- **GitHub**: https://github.com/DerekDardzinski/vaspvis
- **Documentation**: https://derekdardzinski.github.io/vaspvis/
- **PyPI**: https://pypi.org/project/vaspvis/
- **License**: MIT License

## Overview
vaspvis is a highly flexible and customizable Python library for visualizing electronic structure data from VASP calculations. It provides an intuitive API for generating publication-quality band structure and density of states plots with extensive customization options.

**Scientific domain**: Electronic structure visualization, VASP post-processing
**Target user community**: VASP users requiring flexible band structure and DOS visualization

## Theoretical Background
vaspvis processes VASP output to visualize:
- Kohn-Sham eigenvalues along k-paths
- Orbital-projected band character (fatbands)
- Total and projected density of states
- Spin-polarized electronic structure

## Capabilities (CRITICAL)
- **Band Structure**: Standard, projected, and spin-polarized plots
- **DOS/PDOS**: Total and atom/orbital-projected density of states
- **Fatbands**: Orbital character visualization with variable line width
- **Spin Textures**: Spin up/down channel visualization
- **Customization**: Colors, line styles, energy ranges, projections
- **Multi-plot**: Combine band structure and DOS in single figure

## Key Strengths

### Flexibility:
- Extensive customization options
- Multiple projection schemes
- Configurable color maps
- Adjustable plot aesthetics

### VASP Integration:
- Direct VASP file parsing
- EIGENVAL, PROCAR support
- KPOINTS path handling
- POSCAR structure reading

### Publication Quality:
- Matplotlib-based rendering
- Vector format export (PDF, SVG)
- Customizable fonts and sizes
- Professional appearance

## Inputs & Outputs
- **Input formats**:
  - EIGENVAL (eigenvalues)
  - PROCAR (projections)
  - KPOINTS (k-path)
  - POSCAR (structure)
  - INCAR (calculation parameters)
  
- **Output data types**:
  - Matplotlib figures
  - PNG, PDF, SVG images
  - Customizable plots

## Interfaces & Ecosystem
- **Python integration**:
  - Matplotlib for visualization
  - NumPy for data handling
  - Object-oriented API

## Installation
```bash
pip install vaspvis
```

## Usage Examples
```python
from vaspvis import Band, Dos

# Plot band structure
band = Band(folder='path/to/vasp/calc')
band.plot_bands()

# Plot with orbital projections
band.plot_atom_orbitals(atoms=[0], orbitals=[0,1,2])

# Plot DOS
dos = Dos(folder='path/to/vasp/calc')
dos.plot_dos()
```

## Performance Characteristics
- **Speed**: Fast file parsing
- **Memory**: Efficient for standard calculations
- **Visualization**: High-quality matplotlib output

## Limitations & Known Constraints
- **VASP-specific**: Only works with VASP output
- **PROCAR required**: Projections need LORBIT setting
- **k-path**: Requires proper KPOINTS setup

## Comparison with Other Tools
- **vs py4vasp**: vaspvis more customizable, py4vasp official
- **vs sumo**: Both flexible, different API styles
- **vs pyprocar**: Similar capabilities, different interfaces
- **Unique strength**: Highly customizable matplotlib integration

## Application Areas
- Band structure visualization
- DOS analysis
- Orbital character analysis
- Spin-polarized systems
- Publication figure generation

## Best Practices
- Use LORBIT=11 or 12 for projections
- Set appropriate energy range
- Choose meaningful color schemes
- Export as vector formats for publications

## Community and Support
- GitHub issue tracker
- Documentation with examples
- Active development

## Verification & Sources
**Primary sources**:
1. Official documentation: https://derekdardzinski.github.io/vaspvis/
2. GitHub repository: https://github.com/DerekDardzinski/vaspvis

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (GitHub, MIT)
- Developer: Derek Dardzinski
- Active development: Regular updates
