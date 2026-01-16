# Jaguar

## Official Resources
- **Homepage**: https://www.schrodinger.com/products/jaguar
- **Documentation**: Available via Schrödinger Suite
- **Developer**: Schrödinger, Inc.
- **License**: Commercial

## Overview
**Jaguar** is a high-performance *ab initio* electronic structure program known for its unique **pseudospectral** algorithms. Part of the Schrödinger materials science suite, it excels at handling large systems and transition metals with high efficiency. It is heavily used in the pharmaceutical and materials industries for pKa prediction, solvation free energies, and robust geometry optimizations.

**Scientific domain**: Drug Discovery, Catalysis, Materials Design.
**Target user community**: Industrial and academic researchers integrated into the Schrödinger ecosystem.

## Theoretical Methods
- **Basis Set**: Gaussian Type Orbitals (GTOs).
- **Algorithm**: Pseudospectral method (combines analytical and numerical grid integration) for $N^3$ scaling.
- **Functionals**: Extensive DFT functional library (B3LYP-D3, M06, etc.).
- **Solvation**: Poisson-Boltzmann finite element solver (PBF).

## Capabilities
- **Properties**: NMR, IR, VCD, UV-Vis spectra.
- **Reactions**: Transition state search, pKa prediction (highly optimized).
- **Environment**: High-quality solvation models.
- **Analysis**: NBO, spin density, partial charges.

## Key Strengths
- **Pseudospectral efficiency**: Faster than conventional analytical integral codes for medium-large systems.
- **Robustness**: Advanced SCF convergence algorithms (robing, pseudospectral assembly) for difficult transition metals.
- **Integration**: Seamless workflow within Maestro (GUI) and Python API.

## Comparison with Other Codes
- **vs Gaussian**: Jaguar is often faster for DFT due to pseudospectral methods and has superior pKa prediction workflows. Gaussian has a broader range of post-HF methods.
- **vs Orca**: Both are efficient; Jaguar is commercial with strong GUI support, while Orca is free for academia.
## Performance Characteristics
- **Scaling**: Pseudospectral algorithms provide better-than-standard scaling (often $O(N^2)$ to $O(N^3)$) for large DFT calculations.
- **Parallelization**: Highly optimized shared-memory (OpenMP) parallelism; intelligent thread throttling maintains efficiency on multi-core nodes.
- **Limits**: Designed for medium-to-large molecules; extremely large systems (> 2-3 days walltime) are generally considered outside its design scope.

## Limitations & Known Constraints
- **Multi-Reference**: Intentionally lacks multi-reference (CASSCF/MRCI) methods, focusing instead on single-reference DFT/HF for practical throughput.
- **Sensitivity**: Transition state searches on bulky complexes can sometimes yield artificially high barriers due to sensitivity in the initial guess or potential energy surface exploration.

## Best Practices
- **Workflows**: Use the pre-packaged automated workflows (e.g., pKa, Conformational Search) which have optimized default settings/protocols.
- **Approximations**: Always enable Pseudospectral (PS) mode for systems > 50 atoms to gain significant speedups with negligible accuracy loss.

## Community and Support
- **Support**: Commercial support via Schrödinger (very responsive).
- **Seminars**: Regular webinars and training sessions provided by vendor.
## Verification & Sources
**Primary sources**:
1.  **Official Website**: [Schrödinger Jaguar](https://www.schrodinger.com/products/jaguar)
2.  **Literature**: Friesner, R. A., et al. "Jaguar: A high-performance quantum chemistry software program..." *J. Comp. Chem.* (2004).

**Verification status**: ✅ VERIFIED
