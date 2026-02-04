# effmass

## Official Resources
- **Homepage**: https://effmass.readthedocs.io/
- **GitHub**: https://github.com/lucydot/effmass
- **Documentation**: https://effmass.readthedocs.io/
- **PyPI**: https://pypi.org/project/effmass/
- **Publication**: L. Whalley et al., J. Open Source Softw. 3, 797 (2018)
- **License**: MIT License

## Overview
effmass is a Python package for calculating effective masses from ab initio band structure calculations. It provides tools for extracting carrier effective masses using various models including parabolic approximation and Kane's non-parabolic model. The package automatically identifies band extrema and handles multiple valleys.

**Scientific domain**: Semiconductor physics, effective mass calculation, carrier transport
**Target user community**: Researchers studying semiconductors, thermoelectrics, and carrier transport properties

## Theoretical Background
effmass implements effective mass calculations based on:
- Parabolic approximation: m* = ℏ²/(d²E/dk²)
- Kane model for non-parabolic bands: E(1 + E/Eg) = ℏ²k²/2m*
- Optical effective mass: 1/m*_opt = 1/m*_e + 1/m*_h
- DOS effective mass: m*_DOS = (m*_1 · m*_2 · m*_3)^(1/3)
- Conductivity effective mass

## Capabilities (CRITICAL)
- **Effective Mass Calculation**: From DFT band structures
- **Parabolic Fitting**: Standard band curvature analysis
- **Non-parabolic Bands**: Kane model for narrow-gap semiconductors
- **Multi-valley**: Handle degenerate and multiple band extrema
- **Automatic Detection**: Find VBM/CBM automatically
- **Optical Mass**: Combined electron-hole effective mass
- **DOS Mass**: Density of states effective mass
- **Anisotropy**: Direction-dependent effective masses

## Key Strengths

### Multiple Models:
- Parabolic approximation
- Kane non-parabolic model
- Polynomial fitting
- Finite difference methods

### Automatic Analysis:
- Band extrema detection
- Valley identification
- Degeneracy handling
- Energy window selection

### Multi-Code Support:
- VASP (vasprun.xml)
- FHI-aims
- CASTEP
- Extensible to other codes

### Publication Quality:
- Matplotlib visualization
- Data export
- Reproducible analysis

## Inputs & Outputs
- **Input formats**:
  - vasprun.xml (VASP)
  - FHI-aims band output
  - CASTEP band files
  
- **Output data types**:
  - Effective mass values
  - Fitting parameters
  - Band curvature data
  - Matplotlib figures

## Interfaces & Ecosystem
- **Python integration**:
  - NumPy for numerical operations
  - SciPy for fitting
  - Matplotlib for visualization
  - pymatgen for structure handling

## Installation
```bash
pip install effmass
```

## Usage Examples
```python
from effmass import inputs, analysis, outputs

# Load data
settings = inputs.Settings(extrema_search_depth=0.025)
data = inputs.DataVasp(vasprun="vasprun.xml", settings=settings)

# Find segments around band edges
segments = analysis.generate_segments(data)

# Calculate effective masses
for segment in segments:
    print(f"Effective mass: {segment.effective_mass}")
    
# Plot
outputs.plot_segments(data, segments)
```

## Performance Characteristics
- **Speed**: Fast band analysis
- **Accuracy**: Multiple fitting methods
- **Robustness**: Handles complex band structures

## Limitations & Known Constraints
- **k-point density**: Requires dense k-point sampling
- **Band crossing**: May struggle with complex crossings
- **Spin-orbit**: Requires careful handling
- **Indirect gaps**: Need appropriate k-path

## Comparison with Other Tools
- **vs effectivemass (AFLOW)**: effmass has more fitting models
- **vs manual fitting**: effmass is automated, reproducible
- **vs sumo**: effmass specialized for effective mass
- **Unique strength**: Kane model, automatic extrema detection

## Application Areas
- Semiconductor characterization
- Thermoelectric materials
- Solar cell materials
- Transport property prediction
- Band engineering
- Carrier mobility estimation

## Best Practices
- Use dense k-point grids near band edges
- Check fitting quality visually
- Compare parabolic and non-parabolic results
- Consider anisotropy for accurate transport
- Validate against experimental data

## Community and Support
- GitHub issue tracker
- JOSS publication for citation
- Documentation with examples
- Active maintenance

## Verification & Sources
**Primary sources**:
1. Official documentation: https://effmass.readthedocs.io/
2. GitHub repository: https://github.com/lucydot/effmass
3. L. Whalley et al., J. Open Source Softw. 3, 797 (2018)

**Confidence**: VERIFIED - Published in JOSS

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, MIT)
- Developer: Lucy Whalley
- Academic citations: JOSS publication
- Active development: Maintained
