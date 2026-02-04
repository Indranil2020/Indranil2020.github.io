# IFermi

## Official Resources
- **Homepage**: https://fermisurfaces.github.io/IFermi/
- **GitHub**: https://github.com/fermisurfaces/IFermi
- **Documentation**: https://fermisurfaces.github.io/IFermi/
- **PyPI**: https://pypi.org/project/ifermi/
- **Publication**: A. Ganose et al., J. Open Source Softw. 6, 3089 (2021)
- **License**: MIT License

## Overview
IFermi is a Python package for Fermi surface generation, analysis, and visualization from ab initio calculations. It provides tools for extracting Fermi surfaces from DFT band structures, computing Fermi surface properties, and creating publication-quality 2D and 3D visualizations. Developed by the Materials Project team, it integrates seamlessly with pymatgen.

**Scientific domain**: Fermi surface analysis, electronic structure visualization, metal physics
**Target user community**: Researchers studying metallic systems, superconductors, topological materials, and transport properties

## Theoretical Background
IFermi analyzes Fermi surfaces based on:
- Fermi surface: Constant energy surface at E = E_F in k-space
- Isosurface extraction using marching cubes algorithm
- Brillouin zone integration for properties
- Spin texture from spin-polarized calculations
- Fermi velocity: v_F = (1/ℏ)∇_k E(k)

## Capabilities (CRITICAL)
- **Fermi Surface Generation**: Extract isosurfaces from band structures
- **2D Slices**: Fermi surface cross-sections along arbitrary planes
- **3D Visualization**: Interactive 3D Fermi surface plots
- **Spin Textures**: Spin-resolved Fermi surfaces and spin projections
- **Property Calculation**: Fermi velocities, areas, volumes
- **Interpolation**: Band structure interpolation for smooth surfaces
- **Multi-band**: Handle multiple bands crossing Fermi level
- **Export**: Save surfaces in various formats

## Key Strengths

### Visualization Options:
- Mayavi for high-quality 3D rendering
- Plotly for interactive web-based plots
- Matplotlib for 2D slices and publication figures
- Customizable colors and projections

### Analysis Capabilities:
- Fermi surface area calculation
- Average Fermi velocity
- Density of states at Fermi level
- Spin texture analysis
- Band-resolved properties

### Multi-Code Support:
- VASP (vasprun.xml, native)
- Quantum ESPRESSO (via BoltzTraP2)
- BXSF format (universal)
- Extensible to other codes

### Integration:
- pymatgen compatibility
- BoltzTraP2 interface
- Command-line and Python API
- Jupyter notebook friendly

## Inputs & Outputs
- **Input formats**:
  - vasprun.xml (VASP)
  - BXSF files (XCrySDen format)
  - BoltzTraP2 interpolated bands
  - pymatgen BandStructure objects
  
- **Output data types**:
  - 3D Fermi surface meshes
  - 2D slice images
  - Property data (velocities, areas)
  - Interactive HTML plots
  - Image files (PNG, SVG, PDF)

## Interfaces & Ecosystem
- **Python integration**:
  - pymatgen for structure/band handling
  - NumPy/SciPy for numerical operations
  - trimesh for mesh operations
  - matplotlib, plotly, mayavi for visualization
  
- **Framework compatibility**:
  - Materials Project workflows
  - atomate workflows
  - Jupyter notebooks
  - High-throughput screening

## Installation
```bash
pip install ifermi
```

With all visualization backends:
```bash
pip install ifermi[mayavi,plotly]
```

## Usage Examples
Command line:
```bash
# Plot 3D Fermi surface
ifermi plot vasprun.xml

# Plot 2D slice along (001) plane
ifermi plot --slice 0 0 1 vasprun.xml

# Calculate properties
ifermi info vasprun.xml
```

Python API:
```python
from ifermi.surface import FermiSurface
from ifermi.interpolate import FourierInterpolator
from pymatgen.io.vasp import Vasprun

# Load band structure
vr = Vasprun("vasprun.xml")
bs = vr.get_band_structure()

# Interpolate and generate Fermi surface
interpolator = FourierInterpolator(bs)
new_bs, velocities = interpolator.interpolate_bands(return_velocities=True)
fs = FermiSurface.from_band_structure(new_bs, velocities=velocities)

# Plot
fs.plot()
```

## Performance Characteristics
- **Speed**: Fast isosurface extraction
- **Memory**: Efficient mesh handling
- **Interpolation**: Fourier interpolation for smooth surfaces
- **Scalability**: Handles complex multi-sheet Fermi surfaces

## Limitations & Known Constraints
- **Metallic systems only**: Requires bands crossing Fermi level
- **Interpolation quality**: Depends on k-point density
- **Spin-orbit**: Requires careful handling
- **Visualization**: 3D backends may need additional setup
- **Large systems**: Memory for very fine meshes

## Comparison with Other Tools
- **vs FermiSurfer**: IFermi is Python-native, better integration
- **vs XCrySDen**: IFermi has more analysis features
- **vs VESTA**: IFermi specialized for Fermi surfaces
- **Unique strength**: pymatgen integration, spin textures, property calculation

## Application Areas
- Metallic electronic structure
- Superconductivity studies
- Topological material analysis
- Transport property understanding
- Fermi surface nesting analysis
- de Haas-van Alphen effect interpretation
- ARPES comparison

## Best Practices
- Use dense k-point grids for smooth surfaces
- Interpolate bands for better resolution
- Check multiple bands near Fermi level
- Verify with experimental data when available
- Use appropriate energy window for interpolation

## Community and Support
- GitHub issue tracker
- Materials Project community
- JOSS publication for citation
- Active development

## Verification & Sources
**Primary sources**:
1. Official documentation: https://fermisurfaces.github.io/IFermi/
2. GitHub repository: https://github.com/fermisurfaces/IFermi
3. A. Ganose et al., J. Open Source Softw. 6, 3089 (2021)

**Confidence**: VERIFIED - Published in JOSS, Materials Project team

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, MIT)
- Developer: Materials Project team (A. Ganose et al.)
- Academic citations: JOSS publication
- Active development: Regular releases
