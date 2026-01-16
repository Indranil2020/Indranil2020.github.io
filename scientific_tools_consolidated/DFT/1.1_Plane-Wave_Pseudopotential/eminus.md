# eminus

## Official Resources
- Homepage: https://wangenau.github.io/eminus/
- Repository: https://github.com/wangenau/eminus
- Documentation: https://wangenau.github.io/eminus/
- License: Apache License 2.0

## Overview
eminus is a Pythonic electronic structure theory code implementing plane-wave density functional theory (DFT). It is built on top of NumPy and SciPy, prioritizing code readability, simplicity, and extensibility. It is widely used for educational purposes and rapid prototyping of DFT concepts.

**Scientific domain**: Education, algorithm development
**Target user community**: Students, self-learners, researchers prototyping new functionals/methods

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets
- Norm-conserving pseudopotentials (GTH)
- Self-Consistent Field (SCF) methods
- Direct minimization algorithms
- Exchange-correlation functionals (LDA, GGA via Libxc)

## Capabilities
- Ground-state total energy
- Electronic density calculations
- Band structure plotting
- Simple geometry optimization
- Field-dependent calculations
- Multi-state DFT

## Key Strengths

### Pythonic Design:
- Pure Python implementation (with C-backend numerical libraries).
- Modular and easy to read.

### Education:
- Excellent documentation explaining the theory alongside the code.
- "SimpleDFT" prototype included for absolute minimalists.

## Inputs & Outputs
- **Input formats**:
  - Python scripts
  - XYZ structure files
  - GTH pseudopotential files
  
- **Output data types**:
  - VTK/Cube files for visualization
  - Text output (energies)
  - Matplotlib plots

## Interfaces & Ecosystem
- **Libxc**: Interfaces via `pylibxc`.
- **SciPy/NumPy**: Fully integrated.

## Performance Characteristics
- **Speed**: Slower than compiled codes for large systems, but NumPy backend ensures reasonable performance for small systems (up to dozens of atoms).
- **Parallelization**: Relies on threaded libraries (BLAS/LAPACK) and simple multiprocessing options.

## Computational Cost
- **Performance**: Python overhead means it is roughly 10-50x slower than optimized C/Fortran codes for the SCF loop.
- **Bottleneck**: FFT operations and Python loops over grid points (if not vectorized) are the main limits.
- **Memory**: NumPy arrays are efficient, but Python object overhead exists.

## Best Practices

### Learning Strategy:
- **Read the Code**: The source code is written to be read like a textbook.
- **Jupyter**: Run small tutorials in Jupyter notebooks to visualize densities interactively.

### Prototyping:
- **Custom Functionals**: Use eminus to test new Exchange-Correlation functional forms before implementing them in complex C++ codes.

## Community and Support
- **Hosting**: [GitHub](https://github.com/wangenau/eminus).
- **Development**: Active single-developer project with community contributions.
- **Docs**: High-quality Sphinx documentation at [wangenau.github.io/eminus](https://wangenau.github.io/eminus/).

## Verification & Sources
**Primary sources**:
1. Official Website: https://wangenau.github.io/eminus/
2. "eminus: A pythonic electronic structure code" (Software paper)

**Confidence**: VERIFIED - Functional and well-documented.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: Education
- Key Feature: Pythonic Design

## Limitations & Known Constraints
- **System Size**: Not intended for large-scale production runs.
- **Features**: Lacks advanced features like PAW, Hubbard U, or relativistic effects found in major codes.

## Comparison with Other Codes
- **vs PySCF**: PySCF uses Gaussian orbitals; eminus uses plane waves.
- **vs DFTK.jl**: DFTK uses Julia; eminus uses Python. Both target similar "modern/readable" niches.

