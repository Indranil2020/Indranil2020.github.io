# PyGW

## Official Resources
- Homepage: https://github.com/lechifflier/PyGW
- Documentation: https://github.com/lechifflier/PyGW#readme
- Source Repository: https://github.com/lechifflier/PyGW
- License: Open Source

## Overview
PyGW is an electronic structure code for performing G0W0 and GW0 quasiparticle calculations on realistic materials. Implemented as a hybrid Fortran/Python code, it bridges computational efficiency with modern scripting capabilities for GW calculations in condensed matter physics.

**Scientific domain**: Quasiparticle band structures, band gap calculations, electronic excitations  
**Target user community**: Condensed matter physicists studying electronic properties of materials

## Theoretical Methods
- G0W0 approximation
- GW0 (eigenvalue self-consistent)
- Many-body perturbation theory
- Screened Coulomb interaction
- Quasiparticle corrections
- Plane-wave / pseudopotential framework

## Capabilities (CRITICAL)
- G0W0 quasiparticle energies
- GW0 self-consistent eigenvalues
- Band structure calculations
- Band gap predictions
- Quasiparticle corrections to DFT
- Bulk materials
- Semiconductor and insulator systems

**Sources**: Official GitHub repository

## Key Strengths

### Fortran/Python Hybrid:
- Computational efficiency from Fortran
- User-friendly Python interface
- Modern workflow integration
- Scriptable calculations

### GW Implementations:
- Standard G0W0 calculations
- GW0 self-consistency
- Proven methodology
- Materials science focus

### Realistic Materials:
- Production-quality calculations
- Plane-wave accuracy
- Pseudopotential efficiency
- Solid-state applications

## Inputs & Outputs
- **Input formats**:
  - DFT wavefunctions and eigenvalues
  - Pseudopotential files
  - Python configuration
  
- **Output data types**:
  - Quasiparticle energies
  - Band structures
  - Band gaps
  - Self-energy data

## Interfaces & Ecosystem
- **DFT Integration**:
  - Requires DFT input (wavefunctions, eigenvalues)
  - Interfaces with plane-wave DFT codes
  
- **Python scripting**:
  - Python driver scripts
  - Workflow automation
  - Post-processing capabilities

## Performance Characteristics
- **Speed**: Fortran computational core
- **Accuracy**: Standard GW precision
- **System size**: Typical GW scaling
- **Memory**: Plane-wave requirements

## Computational Cost
- **G0W0**: Single-shot calculation
- **GW0**: Multiple iterations for self-consistency
- **Scaling**: O(N^4) typical GW scaling

## Limitations & Known Constraints
- **Documentation**: Limited compared to major codes
- **DFT interface**: Requires specific input format
- **Development**: Smaller community

## Comparison with Other Codes
- **vs BerkeleyGW**: PyGW smaller, BerkeleyGW more features
- **vs Yambo**: Different interface and workflow
- **Unique strength**: Fortran/Python hybrid design

## Application Areas

### Semiconductors:
- Band gap calculations
- Quasiparticle corrections
- Electronic structure

### Insulators:
- Accurate band gaps
- Beyond-DFT corrections
- Materials screening

## Community and Support
- Open-source on GitHub
- Academic development
- Limited but growing documentation

## Verification & Sources
**Primary sources**:
1. Official GitHub: https://github.com/lechifflier/PyGW
2. Active development (2024 commits)

**Confidence**: VERIFIED
- GitHub repository: ACCESSIBLE
- Active development: Yes
- Working implementation: Confirmed
