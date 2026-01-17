# Green-MBPT

## Official Resources
- Homepage: https://green-phys.org/
- Documentation: https://green-phys.org/docs/
- Source Repository: https://github.com/Green-Phys/green-mbpt
- License: MIT License

## Overview
Green-MBPT is a many-body perturbation theory solver within the Green software framework. It provides GW weak-coupling simulations and related many-body calculations for electronic structure, designed for integration with the broader Green computational ecosystem.

**Scientific domain**: Many-body perturbation theory, GW calculations, electronic correlations  
**Target user community**: Researchers using the Green framework for correlated electron calculations

## Theoretical Methods
- GW approximation
- Many-body perturbation theory (MBPT)
- Weak-coupling expansion
- Self-energy calculations
- Green's function methods
- Diagrammatic approaches

## Capabilities (CRITICAL)
- GW weak-coupling simulations
- Self-energy calculations
- Quasiparticle properties
- Green's function evaluation
- Integration with Green framework
- Modular MBPT solvers
- Dielectric response

**Sources**: Official GitHub repository, Green-Phys documentation

## Key Strengths

### Green Framework Integration:
- Part of comprehensive Green ecosystem
- Modular design
- Interoperable components
- Consistent interfaces

### Weak-Coupling GW:
- Perturbative treatment
- Efficient implementation
- Standard GW methodology
- Physical transparency

### Modern Design:
- Active development (2024)
- MIT open source
- Clean codebase
- Good documentation

## Inputs & Outputs
- **Input formats**:
  - Green framework data structures
  - PySCF mean-field input
  - Standard electronic structure data
  
- **Output data types**:
  - Self-energy matrices
  - Quasiparticle energies
  - Green's functions
  - Spectral properties

## Interfaces & Ecosystem
- **Green Framework**:
  - green-integrals
  - green-sc (self-consistency)
  - green-grids
  - green-ac (analytic continuation)
  
- **External DFT**:
  - PySCF integration
  - Standard mean-field methods

## Advanced Features

### Modular Solvers:
- Separable components
- Flexible combinations
- Custom workflows

### Self-Consistency:
- Via green-sc module
- Convergence control
- Multiple schemes

## Performance Characteristics
- **Speed**: Efficient implementation
- **Accuracy**: Standard GW precision
- **System size**: Typical GW scaling
- **Parallelization**: Modern parallel support

## Computational Cost
- **GW**: Polynomial scaling
- **Self-consistency**: Additional iterations
- **Memory**: Green's function storage

## Limitations & Known Constraints
- **Framework dependency**: Best within Green ecosystem
- **Standalone use**: Requires Green infrastructure
- **Documentation**: Growing but evolving

## Comparison with Other Codes
- **vs BerkeleyGW**: Green-MBPT modular, BerkeleyGW standalone
- **vs molgw**: Different framework philosophy
- **Unique strength**: Green ecosystem integration, modular design

## Application Areas

### Electronic Structure:
- Quasiparticle corrections
- Band structure improvements
- Spectral properties

### Correlated Materials:
- Weak correlation regime
- Beyond-DFT corrections
- Material screening

## Community and Support
- Open-source MIT license
- Green-Phys organization
- Active GitHub development
- Documentation available

## Verification & Sources
**Primary sources**:
1. Official GitHub: https://github.com/Green-Phys/green-mbpt
2. Green-Phys website: https://green-phys.org/
3. Methodology publications

**Confidence**: VERIFIED
- GitHub repository: ACCESSIBLE
- Active development: Yes (2024)
- Documentation: Available
