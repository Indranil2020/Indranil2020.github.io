# bandplot

## Official Resources
- PyPI: https://pypi.org/project/bandplot/
- Documentation: Included in package
- License: Open source

## Overview
**bandplot** is a Python package for electron band structure, DOS, and phonon band structure/DOS plotting from VASPKIT or phonopy results. It provides automated plotting with two main commands: `bandplot` for electronic band structures and `phononplot` for phonon dispersions.

**Scientific domain**: Band structure and phonon plotting from VASPKIT/phonopy  
**Target user community**: Researchers using VASPKIT or phonopy who need automated publication-quality plots

## Theoretical Methods
- Band structure plotting from VASPKIT output
- DOS plotting from VASPKIT output
- Phonon band structure from phonopy
- Phonon DOS from phonopy
- matplotlib-based visualization

## Capabilities (CRITICAL)
- Electronic band structure plotting
- Electronic DOS plotting
- Phonon band structure plotting
- Phonon DOS plotting
- VASPKIT output parsing
- phonopy output parsing
- PyPI installable

**Sources**: PyPI, GitHub

## Key Strengths

### VASPKIT Integration:
- Direct VASPKIT output parsing
- Standard VASPKIT workflow
- No manual data extraction
- Consistent with VASPKIT conventions

### Dual Electronic/Phonon:
- Electronic band + DOS
- Phonon band + DOS
- Same tool for both
- Consistent style

### PyPI Installable:
- pip install bandplot
- No manual setup
- Standard Python package
- Easy to use

## Inputs & Outputs
- **Input formats**:
  - VASPKIT band data
  - VASPKIT DOS data
  - phonopy band data
  - phonopy DOS data
  
- **Output data types**:
  - Band structure plots
  - DOS plots
  - Phonon band plots
  - Phonon DOS plots

## Interfaces & Ecosystem
- **VASPKIT**: Primary data source
- **phonopy**: Phonon data source
- **Matplotlib**: Visualization
- **Python**: Core language

## Performance Characteristics
- **Speed**: Fast (plotting)
- **Accuracy**: DFT-level
- **System size**: Any
- **Memory**: Low

## Computational Cost
- **Plotting**: Seconds
- **DFT pre-requisite**: Hours (separate)
- **Typical**: Very efficient

## Limitations & Known Constraints
- **VASPKIT dependent**: Requires VASPKIT-processed data
- **Plotting only**: No analysis features
- **Limited customization**: Standard plots
- **Small community**: Niche tool

## Comparison with Other Codes
- **vs sumo**: bandplot is VASPKIT-specific, sumo is general VASP
- **vs plot4dft**: bandplot is VASPKIT+phonopy, plot4dft is raw VASP/QE
- **vs VASPKIT plotting**: bandplot extends VASPKIT plotting capabilities
- **Unique strength**: PyPI-installable band+DOS+phonon plotting from VASPKIT/phonopy output

## Application Areas

### Publication Figures:
- Band structure plots
- DOS plots
- Phonon dispersion plots
- Combined figures

### VASPKIT Workflow:
- Post-VASPKIT plotting
- Standardized output
- Consistent formatting
- Quick publication figures

## Best Practices

### VASPKIT Setup:
- Use VASPKIT to generate band/DOS data first
- Follow VASPKIT conventions
- Use appropriate k-path
- Check data quality before plotting

### Plotting:
- Use default settings for quick plots
- Customize for publication
- Validate against known systems
- Export in vector format

## Community and Support
- PyPI package
- Open source
- Limited documentation
- Niche community

## Verification & Sources
**Primary sources**:
1. PyPI: https://pypi.org/project/bandplot/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (PyPI)
- Documentation: Included in package
- Specialized strength: PyPI-installable band+DOS+phonon plotting from VASPKIT/phonopy output
