# GreenCheetah

## Official Resources
- Source Repository: https://github.com/StxGuy/GreenCheetah
- Documentation: Included in repository
- License: Open source

## Overview
**GreenCheetah** is a Non-Equilibrium Green's Function (NEGF) approach for quantum transport implemented in Fortran and C++/Armadillo. It provides efficient computation of quantum transport properties in nanoscale devices using the NEGF formalism.

**Scientific domain**: NEGF quantum transport, device simulation  
**Target user community**: Researchers simulating quantum transport in nanoscale devices using NEGF

## Theoretical Methods
- Non-equilibrium Green's function (NEGF)
- Tight-binding Hamiltonians
- Self-energy calculation
- Landauer-Büttiker formalism
- Recursive Green's function algorithm
- Transmission function calculation

## Capabilities (CRITICAL)
- NEGF quantum transport calculation
- Transmission function
- Conductance calculation
- Tight-binding model transport
- Recursive Green's function algorithm
- Fortran and C++ implementations
- Efficient sparse matrix operations

**Sources**: GitHub repository

## Key Strengths

### Performance:
- Fortran/C++ implementation
- Armadillo linear algebra
- Sparse matrix optimization
- Recursive algorithm efficiency

### NEGF Framework:
- Full Green's function calculation
- Self-energy for leads
- Open boundary conditions
- Bias-dependent transport

## Inputs & Outputs
- **Input formats**:
  - Hamiltonian matrix files
  - Lead parameters
  - Transport configuration
  
- **Output data types**:
  - Transmission vs energy
  - Conductance
  - Local density of states
  - Current density

## Interfaces & Ecosystem
- **Fortran/C++**: Core computation
- **Armadillo**: Linear algebra library
- **Python**: Possible scripting interface

## Performance Characteristics
- **Speed**: Fast (compiled code)
- **Accuracy**: Depends on Hamiltonian
- **System size**: Thousands of orbitals
- **Memory**: Moderate

## Computational Cost
- **Transmission**: Seconds to minutes
- **Typical**: Efficient

## Limitations & Known Constraints
- **Tight-binding only**: No DFT integration
- **Limited documentation**: Research code
- **Small community**: Research group code
- **No GUI**: Command-line only

## Comparison with Other Codes
- **vs Gollum**: GreenCheetah is Fortran/C++, Gollum is more feature-rich
- **vs Jiezi**: GreenCheetah is compiled, Jiezi is Python with self-consistent Poisson
- **vs Nanodcal**: GreenCheetah is open source TB, Nanodcal is commercial LCAO
- **Unique strength**: Efficient Fortran/C++ NEGF implementation with Armadillo, recursive Green's function

## Application Areas

### Nanoscale Transport:
- Nanowire conductance
- Molecular junctions
- Quantum dot transport
- 2D material devices

### Method Development:
- NEGF algorithm testing
- Recursive GF benchmarks
- Sparse matrix optimization
- Transport method comparison

## Best Practices

### Hamiltonian Setup:
- Use well-parameterized TB models
- Include sufficient lead layers
- Check convergence of self-energies
- Validate against analytical models

### Performance:
- Use sparse matrix mode
- Optimize memory layout
- Profile hot paths
- Compare Fortran vs C++ paths

## Community and Support
- Open source on GitHub
- Research code
- Limited documentation
- Example calculations provided

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/StxGuy/GreenCheetah

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Research code
- Specialized strength: Efficient Fortran/C++ NEGF implementation with Armadillo, recursive Green's function
