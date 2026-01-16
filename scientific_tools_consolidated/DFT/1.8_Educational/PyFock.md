# PyFock

## Official Resources
- **Homepage**: https://github.com/manassharma07/PyFock
- **Source Repository**: https://github.com/manassharma07/PyFock
- **Developer**: Manas Sharma
- **License**: MIT License

## Overview
**PyFock** is a modern, pure-Python implementation of Density Functional Theory (DFT) and Hartree-Fock (HF) methods. Designed with transparency and hackability in mind, it leverages **JIT compilation (Numba)** and optionally **GPU acceleration (CuPy)** to achieve performance comparable to compiled C++ codes while maintaining the simplicity of a Python codebase.

**Scientific domain**: Algorithm Development, Education, Specialized Python-based QM.
**Target user community**: Developers, students, and researchers prototyping new methods.

## Theoretical Methods
- **Methods**: Restricted and Unrestricted Hartree-Fock (RHF/UHF), DFT (RKS/UKS).
- **Functionals**: LDA, GGA (PBE, BLYP) via LibXC.
- **Integrals**: Obara-Saika scheme implemented in Python/Numba.
- **DIIS**: Direct Inversion in the Iterative Subspace for convergence.

## Capabilities
- **Transparency**: Every step of the SCF cycle is visible in readable Python.
- **Acceleration**: Can switch between NumPy (CPU) and CuPy (GPU) backends.
- **Visualization**: Built-in simple GUI for results (PyFock-GUI).

## Key Strengths
- **Educational Value**: Excellent for understanding how a DFT code is built from scratch.
- **Modern Stack**: Demonstrates high-performance Python computing (JIT/GPU).
- **Open Source**: MIT licensed and actively hackable.

## Comparison with Other Codes
- **vs PySCF**: PySCF is the industry standard for Python QC, but it relies on C modules for heavy lifting. PyFock is *pure* Python (facilitated by JIT).
- **vs Psi4**: Psi4 is C++ with Python bindings.
## Performance Characteristics
- **Precision**: Validated decision better than $10^{-7}$ Hartree compared to PySCF.
- **Scaling**: Near-quadratic $O(N^{2.05})$ for ERIs due to Cauchy-Schwarz screening.
- **GPU Speedup**: Can achieve 14x-20x speedup over CPU (4-core) execution using CuPy backend for large systems.

## Limitations & Known Constraints
- **Features**: Currently lacks analytical gradients (geometry optimization is limited/numerical), periodic boundary conditions, and hybrid functionals (roadmap items).
- **Scope**: Primarily for single-point energy and algorithm demonstrations, not for production dynamics or solid-state physics.

## Best Practices
- **Backend**: Use `numpy` backend for small debugging cases (traceable), switch to `cupy` (GPU) for exploring larger system performance.
- **Memory**: Density fitting is enabled by default to save memory; disable only for debugging raw integrals.

## Community and Support
- **Support**: Community-driven via GitHub Issues.
- **Tutorials**: Author provides detailed blog posts (BragitOff.com) explaining the code structure.
## Verification & Sources
**Primary sources**:
1.  **Repository**: [PyFock GitHub](https://github.com/manassharma07/PyFock)
2.  **Author Blog**: [BragitOff.com](https://www.bragitoff.com/)

**Verification status**: ✅ VERIFIED
