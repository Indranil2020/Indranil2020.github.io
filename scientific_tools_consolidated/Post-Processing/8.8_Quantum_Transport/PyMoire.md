# PyMoire

## Official Resources
- Source Repository: https://github.com/mahyar-servati/PyMoire
- Documentation: Included in repository
- License: Open source

## Overview
**PyMoire** is a Python package for tight-binding calculation of twisted bilayer graphene and other moiré systems based on mapped Wannier functions. It calculates band structures and density of states for twisted bilayer systems at any commensurate twist angle.

**Scientific domain**: Moiré materials, twisted bilayer transport, Wannier-based TB  
**Target user community**: Researchers studying electronic properties of twisted bilayer and moiré superlattice systems

## Theoretical Methods
- Tight-binding model for twisted bilayers
- Mapped Wannier functions
- Commensurate twist angle construction
- Moiré superlattice Hamiltonian
- Interlayer coupling
- Band structure of moiré systems

## Capabilities (CRITICAL)
- Band structure of twisted bilayer graphene
- DOS calculation for moiré systems
- Any commensurate twist angle
- Wannier-based hopping parameters
- Interlayer coupling calculation
- Moiré superlattice construction

**Sources**: GitHub repository

## Key Strengths

### Moiré-Specific:
- Purpose-built for twisted bilayers
- Commensurate angle construction
- Moiré superlattice handling
- Interlayer coupling included

### Wannier-Based:
- Ab initio quality hopping
- Mapped Wannier functions
- Systematic improvement
- DFT-consistent parameters

### Flexible:
- Any commensurate angle
- Multiple bilayer systems
- Customizable parameters
- Extensible framework

## Inputs & Outputs
- **Input formats**:
  - Wannier90 Hamiltonian files
  - Twist angle specification
  - Interlayer coupling parameters
  
- **Output data types**:
  - Band structure of moiré system
  - Density of states
  - Moiré Hamiltonian
  - Flat band characterization

## Interfaces & Ecosystem
- **Wannier90**: Hamiltonian extraction
- **Python**: Core language
- **NumPy**: Numerical computation

## Performance Characteristics
- **Speed**: Moderate (large Hamiltonians)
- **Accuracy**: Wannier-level
- **System size**: Thousands of atoms (moiré cells)
- **Memory**: High for small angles

## Computational Cost
- **Band structure**: Minutes to hours
- **DOS**: Minutes
- **Typical**: Moderate

## Limitations & Known Constraints
- **Commensurate only**: Only commensurate angles
- **Bilayer focus**: Primarily bilayer graphene
- **Wannier dependency**: Requires Wannier90 Hamiltonian
- **Memory**: Large for small twist angles
- **Documentation**: Limited

## Comparison with Other Codes
- **vs Kwant**: PyMoire is moiré-specific, Kwant is general transport
- **vs TBPLaS**: PyMoire is moiré, TBPLaS is general TB
- **vs NanoNet**: PyMoire is Wannier-based, NanoNet is SK-based
- **Unique strength**: Tight-binding calculation of twisted bilayer moiré systems with Wannier functions

## Application Areas

### Twisted Bilayer Graphene:
- Magic angle flat bands
- Superconductivity-related bands
- Correlated insulator states
- Twist angle dependence

### Moiré Materials:
- Moiré superlattice bands
- Flat band engineering
- Interlayer coupling effects
- Twist angle optimization

### Quantum Transport:
- Moiré system conductance
- Flat band transport
- Topological properties
- Valley-dependent transport

## Best Practices

### Wannier Setup:
- Use well-localized Wannier functions
- Include sufficient interlayer hopping
- Validate against DFT bands
- Check Wannier spread convergence

### Twist Angle:
- Start with known commensurate angles
- Test convergence with moiré cell size
- Monitor flat band formation
- Compare with continuum model

## Community and Support
- Open source on GitHub
- Research code
- Limited documentation
- Example calculations provided

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mahyar-servati/PyMoire

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Research code
- Specialized strength: Tight-binding calculation of twisted bilayer moiré systems with Wannier functions
