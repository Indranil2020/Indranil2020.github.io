# TeraChem

## Official Resources
- **Homepage**: https://www.petachem.com/
- **Documentation**: https://www.petachem.com/documentation.html
- **Developer**: PetaChem, LLC
- **License**: Commercial

## Overview
**TeraChem** is a general-purpose quantum chemistry software designed from the ground up for **GPU acceleration**. It was one of the first codes to demonstrate that consumer-grade GPUs could outperform supercomputers for certain quantum chemistry tasks. It provides ultra-fast DFT and *ab initio* molecular dynamics (AIMD) for large molecular systems.

**Scientific domain**: Photochemistry, Biochemistry, AIMD, Nanomaterials.
**Target user community**: Researchers needing high-speed execution for medium-to-large molecules (nanoseconds of AIMD).

## Theoretical Methods
- **Basis Set**: Gaussian Type Orbitals (GTOs).
- **Methodology**: Hartree-Fock and DFT (Kohn-Sham).
- **Algorithms**: GPU-optimized integral evaluation.
- **Excited States**: CIS and TDDFT.

## Capabilities
- **Speed**: Orders of magnitude faster than CPU codes for suitable systems.
- **Dynamics**: Efficient *ab initio* molecular dynamics (AIMD) on a single workstation.
- **Optimization**: Geometry optimization and transition state search.
- **Solvation**: Implicit solvent models (PCM).

## Key Strengths
- **GPU Native**: Not a port; the code structure is designed for SIMD parallelism of GPUs.
- **Throughput**: Enables QM studies on protein-sized systems or long dynamics trajectories previously impossible.
- **Interactive**: Can sometimes be fast enough for "interactive" quantum chemistry.

## Comparison with Other Codes
- **vs Gaussian/Q-Chem**: TeraChem is significantly faster for DFT/HF on GPUs but has a narrower feature set (fewer post-HF methods like CCSD(T)).
- **vs GAMESS (GPU)**: TeraChem is widely considered the most mature GPU-first implementation.
## Performance Characteristics
- **Precision**: Uses a mixed-precision mode (single/double) to maximize consumer GPU throughput.
- **Scaling**: Excellent scaling with system size for suitable molecules, often beating CPU clusters.
- **Hardware**: Designed for NVIDIA GPUs (CUDA).

## Limitations & Known Constraints
- **Periodic Boundary Conditions**: Current PBC implementation operates primarily at the **$\Gamma$-point only**. It is not suitable for systems requiring dense k-point sampling (e.g., metals, small unit cells).
- **Integral Limits**: One-particle basis sets are often limited to lower angular momentum (e.g., d-functions) for GPU acceleration.
- **QM/MM Electrostatics**: Long-range electrostatic interactions in periodic QM/MM can be limited compared to specialized CPU codes.

## Best Practices
- **Hardware**: Use gaming-grade GPUs (e.g., GeForce) for cost-effective performance; high-end Tesla cards are supported but cost/performance ratio is often better on consumer cards for single-precision work.
- **Precision**: Be aware of the mixed-precision thresholds (`dprecision`, `sprecision`) if doing highly sensitive energy comparisons.

## Community and Support
- **Support**: Commercial support via PetaChem.
- **Forum**: TeraChem user forum available.
## Verification & Sources
**Primary sources**:
1.  **Official Website**: [PetaChem](https://www.petachem.com/)
2.  **Literature**: Ufimtsev, I. S., & Martinez, T. J. (2009). "Quantum Chemistry on Graphical Processing Units." *CISE*.

**Verification status**: ✅ VERIFIED
