# kspy-tddft

## Official Resources
- Homepage: https://github.com/pwborthwick/kspy-tddft
- Documentation: https://github.com/pwborthwick/kspy-tddft/blob/master/README.md
- Source Repository: https://github.com/pwborthwick/kspy-tddft
- License: Open Source

## Overview
kspy-tddft is a pure Python implementation of Time-Dependent Density Functional Theory supporting both Linear-Response TDDFT (in the Tamm-Dancoff Approximation) and Real-Time TDDFT. It is designed for educational purposes and small molecular systems, featuring Magnus expansion time propagation and Padé approximant spectral analysis.

**Scientific domain**: Molecular excited states, optical properties, ultrafast dynamics  
**Target user community**: Educators, students, and researchers wanting a transparent, hackable TDDFT implementation in Python

## Theoretical Methods
- Kohn-Sham DFT ground state
- Linear-Response TDDFT (LR-TDDFT)
- Tamm-Dancoff Approximation (TDA)
- Real-Time TDDFT (RT-TDDFT)
- Magnus expansion for time propagation (2nd order)
- Electric dipole response (length gauge)
- External fields: Delta-kick and Gaussian waveforms
- Padé approximants for spectral analysis

## Capabilities
- **Ground-state DFT** calculations
- **LR-TDDFT/TDA** excitation energies
- **RT-TDDFT** time propagation
- **Absorption spectra** via Fourier/Padé analysis
- **Response properties** (polarizabilities)
- **Multiple functionals** support
- **Flexible basis sets** (Gaussian-type)

## Key Strengths

### Educational Design:
- Pure Python implementation
- Transparent, readable code
- Suitable for learning TDDFT methods
- Easy to modify and extend

### Multiple TDDFT Flavors:
- Both LR-TDDFT and RT-TDDFT
- Tamm-Dancoff Approximation
- Full TDDFT (if extended)

### Modern Algorithms:
- Magnus expansion for stable propagation
- Padé approximants for accelerated spectra
- Efficient spectral analysis

### Python Ecosystem:
- NumPy/SciPy integration
- Easy visualization with matplotlib
- Jupyter notebook compatible

## Inputs & Outputs
- **Input formats**:
  - Python scripts/API
  - Molecular geometry (internal format)
  - Basis set specifications
  
- **Output data types**:
  - Excitation energies
  - Oscillator strengths
  - Time-dependent dipole moments
  - Absorption spectra
  - Transition densities

## Interfaces & Ecosystem
- **Pure Python**: NumPy, SciPy dependencies
- **Visualization**: Matplotlib compatible
- **Extensible**: Easy to add new functionals/features

## Performance Characteristics
- **Speed**: Suitable for small molecules (educational)
- **Accuracy**: Depends on basis set and functional
- **System size**: ~10-50 atoms practical
- **Memory**: Limited by Python/NumPy

## Computational Cost
- **Educational Scale**: Designed for seconds-to-minutes calculations on laptops for small molecules.
- **RT-TDDFT**: Cost scales linearly with total propagation time (number of steps).
- **Efficiency**: NumPy operations are vectorized but slower than optimized Fortran/C++.


## Limitations & Known Constraints
- **System size**: Small molecules only (educational focus)
- **Performance**: Pure Python (slower than compiled codes)
- **Basis sets**: Limited library (user can add)
- **Periodic systems**: Not supported (molecular only)
- **Production use**: Designed for learning, not large-scale production

## Comparison with Other Codes
- **vs PySCF**: kspy-tddft smaller, educational; PySCF production-ready
- **vs Psi4**: kspy-tddft pure Python, more transparent; Psi4 faster, broader
- **vs CE-TDDFT**: kspy-tddft molecular/Gaussian; CE-TDDFT periodic/plane-wave
- **Unique strength**: Educational transparency, both LR and RT in one package

## Application Areas
- Teaching TDDFT concepts
- Understanding RT vs LR approaches
- Small molecule spectroscopy
- Algorithm prototyping
- Method development and testing

## Best Practices
- Use for learning and prototyping
- Compare with larger codes for validation
- Small time steps for RT stability
- Converge basis set for accuracy

## Community and Support
- Open-source on GitHub
- Author: P.W. Borthwick
- Educational focus with clear documentation
- Example scripts provided

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/pwborthwick/kspy-tddft
2. README with algorithm descriptions

**Confidence**: VERIFIED
- Repository: ACCESSIBLE (GitHub)
- Code: Complete with examples
- Documentation: Good README with theory

**Verification status**: ✅ VERIFIED
