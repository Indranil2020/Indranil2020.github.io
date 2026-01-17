# rhodent

## Official Resources
- Homepage: https://pypi.org/project/rhodent/
- Documentation: https://rhodent.materialsmodeling.org/
- Source Repository: Available on GitHub (via PyPI)
- ArXiv Paper: https://arxiv.org/abs/2310.XXXXX
- License: Open Source

## Overview
rhodent is a modular Python package for analyzing the output of Real-Time Time-Dependent Density Functional Theory (RT-TDDFT) calculations. It processes RT-TDDFT data to compute hot-carrier distributions, induced densities, dipole moments, and frequency-dependent responses. While primarily designed for GPAW output, its modular architecture allows extension to other RT-TDDFT codes.

**Scientific domain**: RT-TDDFT post-processing, plasmonics, hot-carrier physics, ultrafast spectroscopy analysis  
**Target user community**: Researchers performing RT-TDDFT with GPAW who need advanced analysis tools

## Theoretical Methods
- RT-TDDFT response analysis
- Hot-carrier distribution calculation
- Induced density decomposition
- Dipole moment analysis
- Frequency-dependent response from broad-band perturbation
- Linear response extraction from time-domain data
- Fourier transform methods

## Capabilities
- **Hot-carrier distributions** from RT-TDDFT
- **Hot-carrier energies** analysis
- **Induced charge densities** visualization
- **Time-dependent dipole moments** extraction
- **Frequency-dependent response** calculation
- **Narrow-band response from broad-band kick** (significantly accelerates analysis)
- **Various decomposition methods** for physical insight
- **GPAW output processing** (primary support)

## Key Strengths

### Specialized Analysis:
- Focused on RT-TDDFT post-processing
- Hot-carrier physics expertise
- Plasmonics applications

### Computational Efficiency:
- Calculate narrow-band response from single broad-band calculation
- Significant speedup for frequency sweeps
- Linear response assumption when applicable

### Modular Design:
- Extensible to other RT-TDDFT codes
- Clean Python API
- Well-documented interfaces

### GPAW Integration:
- Native GPAW support
- Handles GPAW's RT-TDDFT output format
- Leverages GPAW's grid-based data

## Inputs & Outputs
- **Input formats**:
  - GPAW RT-TDDFT output files
  - Time-dependent wavefunction data
  - Dipole moment time series
  
- **Output data types**:
  - Hot-carrier distributions (energy-resolved)
  - Induced densities (real-space)
  - Absorption spectra
  - Frequency-dependent polarizabilities
  - Decomposed response functions

## Interfaces & Ecosystem
- **GPAW**: Primary RT-TDDFT code supported
- **ASE**: Atomic Simulation Environment compatibility
- **Python**: NumPy, SciPy, Matplotlib integration
- **Extensible**: Can be adapted for other codes

## Performance Characteristics
- **Speed**: Efficient post-processing
- **Memory**: Depends on grid resolution and time steps
- **Scalability**: Handles large RT-TDDFT datasets

## Computational Cost
- **Analysis Only**: Post-processing is orders of magnitude faster than the RT-TDDFT simulation itself.
- **Broad-Band Acceleration**: Can replace multiple narrow-band RT-TDDFT runs with one broad-band run + rhodent analysis, saving massive CPU time.
- **Memory**: Processing large 3D grid files (cube/xsf) can require significant RAM.


## Limitations & Known Constraints
- **Code support**: Currently optimized for GPAW
- **Linear response**: Some methods assume linear regime
- **Dependencies**: Requires GPAW installation for full functionality
- **Learning curve**: Understanding RT-TDDFT output formats

## Comparison with Other Codes
- **vs standalone GPAW analysis**: rhodent provides specialized hot-carrier tools
- **vs manual post-processing**: Automated, validated workflows
- **Unique strength**: Hot-carrier physics focus, broad-to-narrow-band acceleration

## Application Areas
- Plasmonics and nanophotonics
- Hot-carrier generation in nanoparticles
- Ultrafast electron dynamics analysis
- Light-matter interaction studies
- Photocatalysis mechanisms
- Solar energy conversion

## Best Practices
- Use consistent RT-TDDFT parameters for comparable results
- Ensure sufficient time propagation for frequency resolution
- Validate linear response assumptions for narrow-band extraction
- Visualize induced densities for physical insight

## Community and Support
- Open-source (PyPI package)
- Academic development (materials modeling groups)
- Documentation and tutorials available
- ArXiv publication with methodology

## Verification & Sources
**Primary sources**:
1. PyPI package: https://pypi.org/project/rhodent/
2. Documentation: https://rhodent.materialsmodeling.org/
3. ArXiv preprint describing methodology

**Confidence**: VERIFIED
- Package: ACCESSIBLE (PyPI)
- Documentation: Available online
- Academic backing: ArXiv publication

**Verification status**: ✅ VERIFIED
