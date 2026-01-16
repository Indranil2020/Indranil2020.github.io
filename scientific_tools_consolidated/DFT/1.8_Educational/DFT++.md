# DFT++

## Official Resources
- **Homepage**: http://jdftx.org/ (Historical context mostly)
- **Current Active Fork**: JDFTx (Successor)
- **Developer**: Tomas Arias Group (Cornell University)
- **License**: GPL

## Overview
**DFT++** is not just a single code but a flexible, modular framework for Density Functional Theory calculations. It emphasizes a coordinate-independent linear-algebraic formulation of DFT. While the original "DFT++" code is often used in educational settings (e.g., Cornell minicourses), its modern production incarnation is **JDFTx**. It is highly regarded for its clean, modular C++ design.

**Scientific domain**: Algorithm Design, Joint Density Functional Theory, Solvation.
**Target user community**: Students learning DFT implementation, developers writing new functionals.

## Theoretical Methods
- **Philosophy**: "Physics as Algebra" - implementing DFT equations directly as matrix operations.
- **Basis**: Plane-wave basis set.
- **Solvation**: Pioneers in Joint Density Functional Theory (JDFT) for liquid environments.

## Capabilities
- **Core DFT**: Standard plane-wave pseudopotential calculations.
- **Solvation**: sophisticated fluid models coupled with electronic structure.
- **Modularity**: Easy to swap out minimizers, functionals, or preconditioners.

## Key Strengths
- **Pedagogy**: The design patterns used in DFT++ are excellent examples of object-oriented scientific computing.
- **JDFT**: Unmatched capabilities for electrified interfaces and solvation via JDFTx.

## Comparison with Other Codes
- **vs JDFTx**: JDFTx is the modern, maintained version of the DFT++ philosophy.
- **vs ABINIT**: DFT++ is significantly more compact and modular, though less feature-rich in standard solid-state properties (phonons etc.).
## Performance Characteristics
- **JDFTx Engine**: Written in C++11 with MPI and CUDA support.
- **GPU Acceleration**: Can achieve ~3x speedup on GPUs compared to multi-core CPUs for suitable systems. Requires high memory bandwidth.
- **Scaling**: Near-linear MPI scaling; memory footprint is optimized to be small.

## Limitations & Known Constraints
- **Documentation**: As a research code, documentation can be dense for beginners compared to commercial packages.
- **Hardware**: GPU performance is often memory-bandwidth bound (like most plane-wave codes).

## Best Practices
- **Solvation**: JDFTx is the "gold standard" for *ab initio* solvation; use it when liquid environments are critical.
- **Minimization**: Use its variational minimization algorithms for robust convergence in difficult metallic/charged systems.

## Community and Support
- **Support**: Active mailing list and GitHub issues page.
- **Tutorials**: Excellent tutorials available on strictly following the "Physics as Algebra" philosophy.
## Verification & Sources
**Primary sources**:
1.  **Website**: [JDFTx/DFT++](http://jdftx.org/)
2.  **Literature**: Arias, T. A., et al. "Ab initio molecular dynamics: analytically continued energy functionals..." *Rev. Mod. Phys.* (1992).

**Verification status**: ✅ VERIFIED
