# erlabpy

## Official Resources
- Source Repository: https://github.com/kmnhan/erlabpy
- Documentation: https://erlabpy.readthedocs.io/
- License: Open source

## Overview
**erlabpy** is a complete Python workflow for angle-resolved photoemission spectroscopy (ARPES) experiments. It provides tools to handle, manipulate, and visualize data from ARPES experiments with comprehensive analysis capabilities including momentum conversion, Fermi surface mapping, and band structure analysis.

**Scientific domain**: ARPES workflow, data analysis, visualization  
**Target user community**: Researchers performing ARPES experiments needing a complete analysis workflow

## Theoretical Methods
- ARPES data handling and manipulation
- Momentum conversion
- Fermi surface mapping
- Band dispersion analysis
- Symmetrization and averaging
- Background subtraction
- Self-energy analysis

## Capabilities (CRITICAL)
- Complete ARPES workflow
- Data loading from multiple beamlines
- Momentum conversion
- Fermi surface mapping
- Band dispersion analysis
- Symmetrization
- Self-energy extraction
- Interactive visualization

**Sources**: GitHub repository, ReadTheDocs

## Key Strengths

### Complete Workflow:
- End-to-end ARPES analysis
- From raw data to publication
- Multiple beamline support
- Consistent workflow

### Advanced Analysis:
- Self-energy extraction
- Symmetrization
- Momentum distribution curves
- Energy distribution curves
- Fermi surface mapping

### Interactive:
- Interactive visualization
- Jupyter integration
- Real-time exploration
- Publication-quality output

## Inputs & Outputs
- **Input formats**:
  - Multiple beamline formats
  - ARPES data files
  - Energy/Angle grids
  
- **Output data types**:
  - k-converted data
  - Fermi surface maps
  - Band dispersions
  - Self-energy plots

## Interfaces & Ecosystem
- **NumPy/Xarray**: Data handling
- **Matplotlib/Holoviews**: Visualization
- **Python**: Core language

## Performance Characteristics
- **Speed**: Fast (data processing)
- **Accuracy**: Experimental resolution
- **System size**: Large ARPES datasets
- **Memory**: Moderate to high

## Computational Cost
- **Analysis**: Seconds to minutes
- **No DFT needed**: Experimental data
- **Typical**: Efficient

## Limitations & Known Constraints
- **Experimental data only**: Not for DFT simulation
- **Specific beamline formats**: May need adaptation
- **Memory**: Large datasets can be demanding
- **Learning curve**: Comprehensive tool

## Comparison with Other Codes
- **vs PyARPES**: erlabpy is more complete workflow
- **vs peaks**: erlabpy is comprehensive, peaks is modern framework
- **vs arpespythontools**: erlabpy is full workflow, arpespythontools is lightweight
- **Unique strength**: Complete ARPES workflow with self-energy analysis and interactive visualization

## Application Areas

### ARPES Experiments:
- Full data analysis pipeline
- Multiple beamline support
- Automated processing
- Publication preparation

### Strongly Correlated Systems:
- Self-energy analysis
- Quasiparticle dispersion
- Spectral weight transfer
- Kink analysis

### Topological Materials:
- Dirac cone mapping
- Surface state identification
- Fermi surface topology
- Spin-resolved ARPES

## Best Practices

### Data Loading:
- Use appropriate beamline loader
- Check energy/angle calibration
- Verify Fermi level
- Apply momentum conversion

### Analysis:
- Use symmetrization for gap measurement
- Extract self-energy carefully
- Compare with DFT for validation
- Use interactive tools for exploration

## Community and Support
- Open source on GitHub
- ReadTheDocs documentation
- Active development
- Research community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/kmnhan/erlabpy
2. Documentation: https://erlabpy.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- Specialized strength: Complete ARPES workflow with self-energy analysis and interactive visualization
