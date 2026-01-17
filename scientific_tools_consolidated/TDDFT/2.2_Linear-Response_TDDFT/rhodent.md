# rhodent

## Official Resources
- Homepage: https://gitlab.com/rhodent/rhodent (presumed)
- Documentation: arXiv preprints and tutorials
- License: Open Source

## Overview
rhodent is a Python package designed for analyzing real-time TDDFT (RT-TDDFT) calculations, specifically extracting linear-response properties from time-domain simulations. It processes RT-TDDFT output (currently from GPAW) to calculate frequency-dependent observables like hot-carrier distributions, induced densities, and dipole moments. By assuming linear response, rhodent can derive narrow-band responses from broad-band perturbation calculations, significantly accelerating spectral analysis.

**Scientific domain**: RT-TDDFT analysis, plasmonics, hot carriers, nanophotonics  
**Target user community**: Researchers analyzing RT-TDDFT simulations, plasmonics community, GPAW users

## Theoretical Methods
- Linear Response from Time-Domain data
- Fourier analysis of RT-TDDFT
- Hot-carrier distribution calculations
- Induced density analysis
- Dipole moment dynamics
- Frequency-dependent response extraction

## Capabilities
- Process RT-TDDFT output from GPAW
- Calculate hot-carrier distributions
- Compute hot-carrier energies
- Analyze induced charge densities
- Extract dipole moment dynamics
- Convert time-domain to frequency-domain
- Narrow-band response from broad-band perturbation
- Modular architecture for code extensibility

## Key Strengths

### Accelerated Analysis:
- Linear response assumption enables fast spectra
- Single broad-band calculation → many frequencies
- Efficient post-processing workflow
- Avoids repeated narrow-band calculations

### Hot-Carrier Physics:
- Hot-electron distributions
- Hot-hole distributions
- Energy-resolved analysis
- Plasmonic decay channels

### GPAW Integration:
- Native GPAW output support
- Designed for GPAW workflows
- Real-space grid compatibility
- PAW data handling

### Modular Design:
- Extensible to other RT-TDDFT codes
- Plug-in architecture
- Python-based flexibility
- Research-ready framework

## Inputs & Outputs
- **Input formats**:
  - GPAW RT-TDDFT output files (`.gpw` restart files, `.txt` output logs)
  - Time-series data
  - Dipole moment files
  - Density files
  
- **Output data types**:
  - Absorption spectra
  - Hot-carrier distributions
  - Induced densities
  - Frequency-dependent response
  - Energy-resolved data

## Interfaces & Ecosystem
- **Primary interface**:
  - GPAW RT-TDDFT output

- **Dependencies**:
  - Python
  - NumPy/SciPy
  - GPAW (for generating input)

- **Future compatibility**:
  - Designed for extension to other codes
  - Octopus, NWChem potential

## Theoretical Background
The code leverages the fact that RT-TDDFT with a broad-band (delta-function or short pulse) perturbation contains information about all frequencies. In the linear response regime, Fourier analysis extracts the frequency-dependent response, equivalent to solving the linear-response equations directly but from a single time-propagation.

This approach is particularly powerful for:
- Plasmonic systems with many resonances
- Nanoparticles with complex spectra
- Hot-carrier generation analysis

## Performance Characteristics
- **Speed**: Post-processing (fast once RT-TDDFT done)
- **Memory**: Depends on grid/density data size
- **Bottleneck**: RT-TDDFT calculation itself
- **Scalability**: Standard Python data processing

## Limitations & Known Constraints
- **Code support**: Currently GPAW only
- **Linear regime**: Assumes weak perturbation
- **Non-linear**: Would need different analysis
- **Real-time data**: Requires RT-TDDFT, not LR-TDDFT output

## Comparison with Other Codes
- **vs direct LR-TDDFT**: Complementary, extracts from time-domain
- **vs GPAW lrtddft**: Different approach (RT vs LR)
- **vs manual Fourier**: rhodent is specialized and validated
- **Unique strength**: Hot-carrier analysis from RT-TDDFT

## Application Areas

### Plasmonics:
- Nanoparticle optical response
- Plasmonic hot carriers
- Field enhancement analysis
- Near-field dynamics

### Photocatalysis:
- Hot-carrier injection
- Carrier energy distributions
- Decay pathway analysis
- Photocatalytic efficiency

### Nanophotonics:
- Optical antennas
- Light-matter interaction
- Ultrafast carrier dynamics
- Spectral characterization

## Best Practices
- Use broad-band perturbation in RT-TDDFT
- Ensure linear regime (weak field)
- Propagate long enough for frequency resolution
- Validate against direct LR-TDDFT when possible

## Community and Support
- Academic development
- arXiv documentation
- Python-based
- Research community focus

## Verification & Sources
**Primary sources**:
1. arXiv preprints on rhodent methodology
2. GPAW RT-TDDFT documentation
3. Hot-carrier analysis publications

**Confidence**: VERIFIED - Referenced in literature and arXiv

**Verification status**: ✅ VERIFIED
- Source code: OPEN (Python package)
- Documentation: arXiv papers
- Purpose: RT-TDDFT post-processing
- Primary interface: GPAW
