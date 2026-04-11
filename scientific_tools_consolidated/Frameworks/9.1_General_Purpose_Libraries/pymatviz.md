# pymatviz

## Official Resources
- Source Repository: https://github.com/janosh/pymatviz
- Documentation: https://pymatviz.readthedocs.io/
- PyPI: https://pypi.org/project/pymatviz/
- License: Open source (MIT)

## Overview
**pymatviz** is a Python toolkit for visualizations in materials informatics. It provides publication-quality plotting functions for common materials science data types including parity plots, histogram plots, crystal structure visualizations, and phase diagrams, built on matplotlib and plotly.

**Scientific domain**: Materials informatics visualization, ML model analysis, phase diagrams  
**Target user community**: Researchers needing publication-quality plots for materials ML and property analysis

## Theoretical Methods
- Parity and residual plots for ML model evaluation
- Histogram and cumulative distribution plots
- Crystal structure visualization
- Phase diagram plotting
- Sunburst charts for element distributions
- Spacegroup histograms
- Phonon band structure plotting

## Capabilities (CRITICAL)
- Parity plots (with R², MAE, RMSE annotations)
- Residual and histogram plots
- Crystal structure 2D/3D visualization
- Phase diagram plotting
- Spacegroup distribution charts
- Element frequency sunburst charts
- Phonon band structure plotting
- Matplotlib and Plotly backends

**Sources**: GitHub repository, ReadTheDocs

## Key Strengths

### Publication Quality:
- Consistent styling
- Automatic annotations (R², MAE, RMSE)
- Matplotlib and Plotly support
- Customizable themes

### Materials-Specific:
- Crystal structure plots
- Phase diagrams
- Spacegroup histograms
- Element distribution charts
- Phonon band structures

### Easy to Use:
- PyPI installable
- One-function plotting
- Sensible defaults
- Flexible customization

## Inputs & Outputs
- **Input formats**:
  - NumPy arrays
  - pandas DataFrames
  - pymatgen structures
  - phonopy data
  
- **Output data types**:
  - Publication-quality plots
  - SVG/PNG/PDF export
  - Interactive Plotly figures

## Interfaces & Ecosystem
- **pymatgen**: Structure handling
- **matplotlib**: Static plotting
- **plotly**: Interactive plotting
- **pandas**: Data handling

## Performance Characteristics
- **Speed**: Fast (plotting)
- **Accuracy**: N/A (visualization)
- **System size**: Any
- **Memory**: Low

## Computational Cost
- **Plotting**: Seconds
- **No DFT needed**: Post-processing only
- **Typical**: Very efficient

## Limitations & Known Constraints
- **Visualization only**: No analysis features
- **Python 3.9+**: Recent Python required
- **Plotly dependency**: Optional but large
- **Niche**: Materials science specific

## Comparison with Other Codes
- **vs matplotlib directly**: pymatviz is materials-specific, auto-annotated
- **vs sumo**: pymatviz is ML/property focused, sumo is band structure
- **vs matminer**: pymatviz is visualization, matminer is data mining
- **Unique strength**: Materials-specific publication-quality visualization with automatic ML metrics annotation

## Application Areas

### ML Model Evaluation:
- Parity plots for regression
- Residual analysis
- Model comparison
- Feature importance

### Materials Analysis:
- Phase diagram visualization
- Spacegroup statistics
- Element distribution
- Crystal structure display

### Publication:
- Consistent figure styling
- Auto-annotated plots
- High-quality export
- Interactive and static options

## Best Practices

### Plotting:
- Use consistent styling across paper
- Include R², MAE, RMSE annotations
- Export in vector format (SVG/PDF)
- Use interactive Plotly for exploration

### Data:
- Clean data before plotting
- Use appropriate plot types
- Validate against known systems
- Label axes clearly

## Community and Support
- Open source (MIT)
- PyPI installable
- ReadTheDocs documentation
- Active development
- GitHub: janosh/pymatviz

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/janosh/pymatviz
2. Documentation: https://pymatviz.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- PyPI: AVAILABLE
- Specialized strength: Materials-specific publication-quality visualization with automatic ML metrics annotation
