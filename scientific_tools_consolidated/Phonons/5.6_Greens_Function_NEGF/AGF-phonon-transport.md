# AGF-phonon-transport

## Official Resources
- Homepage: https://github.com/brucefan1983/AGF-phonon-transport
- Source Repository: https://github.com/brucefan1983/AGF-phonon-transport
- License: Open Source

## Overview
AGF-phonon-transport is a MATLAB implementation of the Atomistic Green's Function (AGF) method for phonon transport calculations. It enables computation of phonon transmission and thermal conductance across interfaces and nanostructures.

**Scientific domain**: Phonon transport, Green's function methods, interfacial thermal conductance  
**Target user community**: Researchers studying phonon transport at interfaces

## Theoretical Methods
- Atomistic Green's Function (AGF)
- Phonon transmission function
- Landauer formalism
- Surface Green's functions
- Interfacial thermal conductance
- Ballistic phonon transport

## Capabilities (CRITICAL)
- Phonon transmission calculation
- Interfacial thermal conductance
- Green's function computation
- Surface self-energies
- Frequency-resolved transport
- MATLAB implementation

## Key Strengths

### AGF Method:
- Exact for ballistic transport
- Interface-focused
- Frequency resolution
- Quantum effects

### MATLAB Implementation:
- Easy to understand
- Educational value
- Modifiable code
- Quick prototyping

## Inputs & Outputs
- **Input formats**:
  - Force constants
  - Structure information
  - Interface geometry
  
- **Output data types**:
  - Transmission function
  - Thermal conductance
  - Green's functions

## Interfaces & Ecosystem
- **MATLAB**: Primary platform
- Standalone code


## Advanced Features
- **Atomistic Green's Function**: Exact ballistic phonon transport
- **Surface self-energies**: Lead-device coupling
- **Transmission function**: Frequency-resolved transport
- **Landauer formalism**: Thermal conductance calculation
- **Interface modeling**: Heterojunction phonon transport
- **Educational code**: Well-commented MATLAB implementation

## Performance Characteristics
- MATLAB-based: Moderate speed
- Matrix operations: Efficient for small systems
- Memory: Scales with system size

## Computational Cost
- Force constants: External (DFT or empirical)
- AGF calculation: Minutes to hours
- Scales with interface size and frequency points
- Overall: Moderate for typical interface systems

## Best Practices
- Validate force constants before transport calculation
- Check convergence with lead size
- Use sufficient frequency resolution
- Compare with analytical limits when available

## Limitations & Known Constraints
- MATLAB required
- Ballistic limit only
- Limited system sizes
- Educational focus

## Application Areas
- Interfacial thermal resistance
- Nanostructure transport
- Phonon engineering
- Thermal interface materials

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/brucefan1983/AGF-phonon-transport

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
