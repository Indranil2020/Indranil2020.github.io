# tinydft

## Official Resources
- Homepage: https://github.com/theochem/tinydft
- Source Repository: https://github.com/theochem/tinydft
- License: MIT License

## Overview
tinydft is a minimalistic implementations of Density Functional Theory (DFT) in Python, designed specifically for educational purposes. It illustrates the core components of a DFT calculation—grids, basis sets, Hamiltonian construction, and self-consistency—within a compact codebase (often fitting in a single file or small module). It serves as a pedagogical tool for students and developers to read and understand the "black box" of quantum chemistry software.

**Scientific domain**: Educational Quantum Chemistry, Algorithmic DFT
**Target user community**: Students, Teachers, Python developers

## Theoretical Methods
- **Kohn-Sham DFT**: Standard formulation.
- **Exchange-Correlation**: LDA (Local Density Approximation) typically implemented via PySCF or LibXC interfaces, or simple hardcoded Slater exchange.
- **Basis Sets**: Gaussian Type Orbitals (GTOs) primarily.
- **Numerical Integration**: Becke grids for XC integration.
- **DIIS**: Direct Inversion in the Iterative Subspace for convergence acceleration.

## Capabilities (CRITICAL)
- **SP Calculation**: Single Point energy calculation.
- **Density Visualization**: Evaluation of electron density on grids.
- **Matrices**: Exposure of Core Hamiltonian, J (Coulomb), and K (Exchange) matrices.
- **Modifiable Loop**: Easy customization of the SCF cycle.

## Key Strengths

### Simplicity:
- Focuses on readability over performance. 
- Uses standard NumPy syntax to express linear algebra directly from textbook equations.

### Dependency Minimalist:
- Usually requires only NumPy/SciPy (and sometimes PySCF for integrals).
- Easy to install and run anywhere.

## Inputs & Outputs
- **Inputs**:
  - Molecular geometry (XYZ or internal).
  - Basis set name (e.g., 'sto-3g', '6-31g').
- **Outputs**:
  - Total Energy.
  - Orbital Energies (HOMO/LUMO).
  - Density Matrix.

## Interfaces & Ecosystem
- **Python**: Pure Python.
- **PySCF**: Often uses PySCF's `gto` module to handle the complex task of Gaussian integral evaluation, focusing the code on the DFT logic itself.

## Advanced Features
- **XC Functionals**: Can demonstrate implementation of custom functionals by modifying the integration kernel.

## Performance Characteristics
- **Speed**: Slow. Python loops over grids are not optimized for production.
- **Scale**: Limited to very small molecules (H2O, CH4) for demonstration.

## Computational Cost
- **Educational**: Negligible for intended small systems; exponential scaling if pushed to large systems due to unoptimized routines.

## Limitations & Known Constraints
- **Not for Research**: Do not use for publication-quality data.
- **Features**: Lacks geometry optimization, frequencies, or advanced properties.

## Comparison with Other Codes
- **vs PySCF**: PySCF is a production code with C backends; tinydft is a toy code for *learning* how PySCF works.
- **vs Psi4NumPy**: Similar goals; Psi4NumPy utilizes the C++ Psi4 core; tinydft often tries to be more self-contained or relies only on PySCF integrals.
- **Unique strength**: Extreme minimalism for identifying the atomic components of a DFT code.

## Application Areas
- **Classroom**: Interactive coding of a DFT SCF loop in a 1-hour workshop.
- **Algorithm Design**: Testing a new mixing scheme in 50 lines of code.

## Best Practices
- **Read the Source**: The value is in reading the code, not just running it.
- **Use Small Basis**: Stick to STO-3G or 3-21G to keep fast execution.

## Community and Support
- **GitHub**: Educational repo maintenance.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/theochem/tinydft

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Purpose: purely educational.
