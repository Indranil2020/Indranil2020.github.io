# emc

## Official Resources
- Source Repository: https://github.com/afonari/emc
- Documentation: Included in repository
- License: Open source

## Overview
**emc** (Effective Mass Calculator) is a program for calculating effective masses at band extrema in semiconductors using the finite difference method. Available in both FORTRAN and Python versions, it works with VASP and Quantum ESPRESSO outputs to compute anisotropic effective mass tensors.

**Scientific domain**: Effective mass calculation, semiconductor band structure  
**Target user community**: Researchers computing carrier effective masses from DFT band structures

## Theoretical Methods
- Finite difference method for effective masses
- Second derivative of E(k) at band extrema
- Inverse effective mass tensor
- Anisotropic effective mass calculation
- VASP and QE band structure input

## Capabilities (CRITICAL)
- Effective mass calculation at band extrema
- Anisotropic effective mass tensor
- Finite difference approach
- VASP output parsing
- Quantum ESPRESSO output parsing
- FORTRAN and Python implementations
- Principal effective masses

**Sources**: GitHub repository, Comput. Phys. Commun.

## Key Strengths

### Multi-Code Support:
- VASP output parsing
- Quantum ESPRESSO output parsing
- FORTRAN version (fast)
- Python version (flexible)

### Anisotropic Masses:
- Full effective mass tensor
- Principal components
- Direction-dependent masses
- Highly anisotropic materials

### Finite Difference:
- Simple and robust
- No perturbation theory needed
- Works with any DFT code output
- Well-tested

## Inputs & Outputs
- **Input formats**:
  - VASP EIGENVAL
  - QE band structure output
  - Band extrema specifications
  
- **Output data types**:
  - Effective mass tensor
  - Principal effective masses
  - Direction-dependent masses

## Interfaces & Ecosystem
- **VASP**: DFT backend
- **Quantum ESPRESSO**: DFT backend
- **FORTRAN/Python**: Dual implementation

## Performance Characteristics
- **Speed**: Fast (post-processing)
- **Accuracy**: Good (finite difference)
- **System size**: Any
- **Memory**: Low

## Computational Cost
- **Effective mass**: Seconds
- **DFT pre-requisite**: Hours (separate)
- **Typical**: Very efficient

## Limitations & Known Constraints
- **Finite difference**: Numerical errors possible
- **Band extrema required**: Must identify extrema
- **Parabolic assumption**: Limited to parabolic regions
- **Step size sensitivity**: Results depend on k-spacing

## Comparison with Other Codes
- **vs mstar**: emc uses finite differences, mstar uses perturbation theory (WIEN2k)
- **vs effmass**: emc is FORTRAN/Python, effmass is Python-only (VASP)
- **vs Effective-mass-fitting**: emc is finite difference, fitting is polynomial
- **Unique strength**: Dual FORTRAN/Python, VASP+QE support, anisotropic effective mass tensor

## Application Areas

### Semiconductor Physics:
- Carrier effective mass for transport
- Anisotropic mass for device modeling
- Mass tensor for mobility calculation
- Band structure characterization

### Thermoelectric Materials:
- DOS effective mass
- Conductivity effective mass
- Mass anisotropy effects
- Transport property estimation

### Optoelectronics:
- Reduced mass for excitons
- Effective mass for absorption
- Carrier mobility estimation
- Band gap engineering

## Best Practices

### k-Point Selection:
- Use fine k-grid near extrema
- Test convergence with k-spacing
- Use symmetric k-paths
- Validate against known systems

### Finite Difference:
- Choose appropriate step size
- Test step size convergence
- Compare with perturbation theory
- Validate against experiment

## Community and Support
- Open source on GitHub
- Developed by A. Fonari
- Published methodology
- Research code

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/afonari/emc

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Research code
- Specialized strength: Dual FORTRAN/Python effective mass calculator, VASP+QE support, anisotropic tensor
