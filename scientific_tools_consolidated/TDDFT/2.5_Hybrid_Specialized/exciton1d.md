# exciton1d

## Official Resources
- Homepage: https://github.com/nicholashestand/exciton1d
- Source Repository: https://github.com/nicholashestand/exciton1d
- License: MIT License

## Overview
exciton1d is a specialized software package for simulating exciton dynamics and spectroscopy in one-dimensional molecular aggregates. It implements the Holstein Hamiltonian to treat Frenkel excitons coupled to vibrational modes, as well as charge-transfer states. It is designed to calculate absorption spectra, band dispersion, and exciton coherence properties.

**Scientific domain**: Exciton dynamics, molecular aggregates, Frenkel excitons, spectroscopy
**Target user community**: Researchers studying J-aggregates, H-aggregates, and organic semiconductor models

## Theoretical Methods
- Holstein Hamiltonian
- Frenkel Exciton derivation
- Charge-Transfer (CT) states
- Vibrational coupling (Holstein model)
- Heitler-London approximation
- CES approximation
- Two-particle approximation

## Capabilities (CRITICAL)
- Absorption spectra calculation
- Band structure / Dispersion relations
- Density of States (DOS)
- Exciton coherence length analysis
- Inclusion of vibrational disorder
- Inclusion of static disorder
- Time-dependent properties

**Sources**: GitHub repository

## Key Strengths

### Specialized 1D Models:
- Highly optimized for linear chains
- Treats vibronic coupling explicitly
- Analytic and numerical solutions

### Disorder Handling:
- Gaussian static disorder
- Dynamic disorder enabled
- Ensemble averaging

### CT State Mixing:
- Beyond simple Frenkel model
- Important for organic photovoltaics

## Inputs & Outputs
- **Input formats**:
  - Python scripts / Input configurations
  - Interaction parameters (J couplings)
  
- **Output data types**:
  - Spectra (Absorption/Emission)
  - Wavefunctions
  - Eigenvalues
  - Coherence function

## Interfaces & Ecosystem
- **Language**: Python (NumPy/SciPy)
- **Integration**: Can use parameters derived from QC codes

## Advanced Features

### Multi-particle Basis:
- One- and Two-particle basis sets
- Converged vibronic spectra

### Hamiltonian Construction:
- Flexible definition of site energies
- Nearest-neighbor and long-range coupling

## Performance Characteristics
- **Speed**: Very fast (model Hamiltonian)
- **scaling**: N^2 or N^3 depending on basis truncation
- **System size**: Long chains possible (100s of units)

## Computational Cost
- **Low**: Model Hamiltonian diagonalization
- **High**: If basis size grows with vibrations

## Limitations & Known Constraints
- **Model**: Restricted to 1D aggregate models
- **Parameters**: Requires input parameters (not ab initio)
- **Geometry**: Implicit linear geometry

## Comparison with Other Codes
- **vs MCTDH**: exciton1d is specialized/simplified for aggregates
- **vs Ab initio**: exciton1d uses parameterized models, much faster
- **Unique strength**: Dedicated tool for 1D vibronic exciton models

## Application Areas
- **J-aggregates**: Cyanine dyes
- **Organic Semiconductors**: P3HT chains
- **Photosynthesis**: Antenna complex models

## Best Practices
- **Parameterization**: Derive J and E0 from reliable QC
- **Basis Convergence**: Check particle number limit
- **Disorder**: Sufficient ensemble averaging

## Community and Support
- Open-source MIT
- GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/nicholashestand/exciton1d

**Confidence**: VERIFIED - GitHub project

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (MIT)
- Specialized strength: 1D Frenkel-Holstein exciton model
