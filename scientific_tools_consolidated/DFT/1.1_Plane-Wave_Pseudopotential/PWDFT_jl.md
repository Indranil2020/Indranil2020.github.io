# PWDFT

## Official Resources
- Repository: https://github.com/f-fathurrahman/PWDFT.jl
- Documentation: In-repository examples
- License: MIT License

## Overview
PWDFT.jl is a Julia package designed to solve the Kohn-Sham equations using a plane-wave basis set and pseudopotentials. It focuses on providing a clear implementation of standard Plane-Wave DFT methods, making it suitable for educational purposes and for researchers who wish to understand the implementation details of plane-wave DFT codes.

**Scientific domain**: Electronic structure, education
**Target user community**: Students, educators, researchers exploring Julia for DFT

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets
- Norm-conserving pseudopotentials (GTH)
- Direct minimization and SCF methods
- LDA and GGA functionals (Libxc via Libxc.jl)

## Capabilities
- Total energy calculation
- Electron density generation
- Support for standard GTH pseudopotentials
- Simple geometry optimization
- Band structure calculation

## Key Strengths

### Educational Value:
- Concise codebase (~few thousand lines)
- Demonstrates core components of a PW-DFT code clearly
- Written in Julia for readability and performance

### Simplicity:
- Minimal dependencies
- Straightforward installation

## Inputs & Outputs
- **Input formats**:
  - Julia scripts
  - XYZ structure files
  - GTH pseudopotential format
  
- **Output data types**:
  - Text output (energies, eigenvalues)
  - Simple visualization (XSF export likely)

## Interfaces & Ecosystem
- **Julia**: Standard Julia package structure.
- **Libxc**: Uses Libxc for exchange-correlation.

## Computational Cost
- **Scale**: Limited to small systems (tens of atoms) for reasonable runtimes on CPUs.
- **Memory**: Standard Julia array overheads apply.

## Best Practices

### Educational Use:
- **Source Reading**: The code is designed to be read. Users should examine the `src/` directory to understand how constructing the Hamiltonian works.
- **Small Examples**: Stick to simple molecules (H2, H2O) or simple solids (Si) to learn the workflow without waiting for heavy compute.

### Development:
- **Prototyping**: Good for testing simple changes to the SCF loop structure before moving to complex packages like DFTK.

## Community and Support
- **Maintenance**: Primarily a personal project by F. Fathurrahman; updates are sporadic.
- **Support**: Best via GitHub issues or direct academic contact.
- **Resources**: Repository includes examples/tutorials in Jupyter notebooks.

## Performance Characteristics
- **Speed**: Reasonable for small systems; comparable to other prototype codes.
- **Parallelization**: Basic threading support; not primarily designed for massive HPC scaling.

## Limitations & Known Constraints
- **Features**: Much more limited than minimal production codes (like Abinit or QE).
- **Scope**: Primarily an educational/research tool, not for large-scale production.

## Comparison with Other Codes
- **vs DFTK.jl**: DFTK is a broader community project with more modern features (AD, ecosystem integration). PWDFT.jl is a more standalone/personal implementation often used for teaching.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/f-fathurrahman/PWDFT.jl
2. Developer's companion site/notes (f-fathurrahman).

**Confidence**: VERIFIED - Code exists and is functional.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: DFT/Julia
- Key Feature: Educational Implementation
