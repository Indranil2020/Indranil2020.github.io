# momentGW

## Official Resources
- Homepage: https://github.com/BoothGroup/momentGW
- Documentation: https://github.com/BoothGroup/momentGW#readme
- Source Repository: https://github.com/BoothGroup/momentGW
- PyPI: https://pypi.org/project/momentGW/
- License: Apache License 2.0

## Overview
momentGW is a Python package for GW approximation calculations using moment-conserving solutions to the Dyson equation. Built on PySCF, it provides an efficient framework for calculating quasiparticle energies and spectral properties with novel algorithmic advantages including exact frequency integration and avoided analytical continuation.

**Scientific domain**: Quasiparticle energies, ionization potentials, electron affinities, band structures  
**Target user community**: Researchers needing efficient GW calculations with PySCF ecosystem integration

## Theoretical Methods
- GW approximation (G0W0, evGW, qsGW)
- Moment-conserving reformulation of GW theory
- Self-consistent GW schemes
- Dyson equation exact solution
- dTDA and dRPA polarizabilities
- Tensor hypercontraction support
- Resolution of Identity (RI) approximation

## Capabilities (CRITICAL)
- G0W0 quasiparticle calculations
- Eigenvalue self-consistent GW (evGW)
- Quasiparticle self-consistent GW (qsGW)
- Unrestricted (spin-polarized) GW
- Periodic boundary conditions
- Contour deformation methods
- Analytic continuation (avoided via moment approach)
- Spectral function calculations
- Ionization potentials and electron affinities
- Molecular and periodic systems

**Sources**: Official GitHub repository, published methodology papers

## Key Strengths

### Moment-Conserving Approach:
- Efficient starting point calculations
- No approximations in frequency integration
- Exact Dyson equation solution
- Avoids analytical continuation errors
- No iterative quasiparticle equation solving

### PySCF Integration:
- Seamless interface with PySCF
- Access to PySCF mean-field methods
- Compatible with PySCF infrastructure
- Familiar API for PySCF users
- Extensible design

### Self-Consistency Options:
- Multiple self-consistent GW schemes
- Easy switching between methods
- Flexible workflow design
- Consistent implementation

### Modern Implementation:
- Python-based with NumPy/SciPy
- Clean object-oriented design
- Active development (2024)
- Well-documented codebase
- Open-source Apache 2.0

## Inputs & Outputs
- **Input formats**:
  - PySCF molecule/cell objects
  - Basis set specifications
  - Mean-field reference (HF/DFT)
  
- **Output data types**:
  - Quasiparticle energies
  - Self-energy matrices
  - Spectral functions
  - Ionization potentials
  - Electron affinities

## Interfaces & Ecosystem
- **Python integration**:
  - PySCF core dependency
  - NumPy/SciPy compatibility
  - Standard Python ecosystem
  
- **Framework integrations**:
  - PySCF mean-field methods
  - Built-in MP2/RPA interfaces
  - Extendable to other post-HF

## Advanced Features

### Tensor Hypercontraction:
- Reduced memory scaling
- Accelerated tensor operations
- Efficient for larger systems

### Periodic Calculations:
- Full k-point support
- Brillouin zone sampling
- Solid-state applications

### Spin-Polarization:
- Unrestricted reference support
- Open-shell systems
- Magnetic systems

## Performance Characteristics
- **Speed**: Efficient moment-conserving algorithm
- **Accuracy**: Exact frequency integration
- **System size**: Medium molecules to small solids
- **Memory**: Optimized with tensor hypercontraction
- **Parallelization**: NumPy parallel backends

## Computational Cost
- **G0W0**: Efficient single-shot
- **evGW**: Multiple iterations
- **qsGW**: Highest cost, highest accuracy
- **Scaling**: Polynomial with system size

## Limitations & Known Constraints
- **PySCF dependency**: Requires PySCF installation
- **Python overhead**: Some overhead vs Fortran codes
- **Large systems**: Best for small-medium systems
- **Development**: Newer code, evolving features

## Comparison with Other Codes
- **vs molgw**: momentGW uses moment-conserving, molgw traditional approach
- **vs BerkeleyGW**: momentGW molecular focus, BerkeleyGW plane-wave solids
- **vs PySCF-GW**: Different algorithmic approach
- **Unique strength**: Moment-conserving exact frequency integration

## Application Areas

### Molecular Spectroscopy:
- Ionization potentials
- Electron affinities
- HOMO-LUMO gaps
- Photoemission simulation

### Materials Science:
- Small periodic systems
- Band gap predictions
- Electronic structure

## Best Practices

### Reference Selection:
- Test HF vs DFT starting points
- Consider evGW/qsGW for starting point independence
- Verify convergence with basis set

### Convergence:
- Check basis set convergence
- Verify k-point convergence for periodic
- Test self-consistency options

## Community and Support
- Open-source Apache 2.0
- Active GitHub development
- Booth Group (King's College London) maintained
- Documentation and examples available
- Academic publications describing methodology

## Verification & Sources
**Primary sources**:
1. Official GitHub: https://github.com/BoothGroup/momentGW
2. O. J. Backhouse, G. H. Booth, J. Chem. Theory Comput. (methodology papers)
3. PyPI package page

**Confidence**: VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: AVAILABLE
- Active development: Yes (2024)
- Academic papers: Published methodology
