# BAND (Amsterdam Modeling Suite)

## Official Resources
- **Homepage**: https://www.scm.com/product/band/
- **Documentation**: https://www.scm.com/doc/BAND/
- **Developer**: Software for Chemistry & Materials (SCM)
- **License**: Commercial (Academic pricing available)

## Overview
**BAND** is the periodic Density Functional Theory (DFT) code within the Amsterdam Modeling Suite (AMS). Unlike most periodic codes that use plane waves (like VASP or QE), BAND utilizes **atom-centered numerical orbitals (STOs/NAOs)**. This basis set allows for an accurate treatment of both core and valence electrons and makes the code particularly efficient for low-dimensional systems (1D polymers, 2D slabs) and empty space.

**Scientific domain**: Periodic conventions, Surface Science, Catalysis, Nanotubes.
**Target user community**: Chemists and materials scientists dealing with periodic systems where chemical insight (orbitals, bonds) is paramount.

## Theoretical Methods
- **Basis Set**: Slater-Type Orbitals (STOs) and Numerical Atomic Orbitals (NAOs).
- **Functionals**: Extensive library from LDA/GGA to hybrid functionals and dispersion corrections.
- **Relativity**: Scalar relativistic (ZORA) and spin-orbit coupling.
- **Solvation**: COSMO model for periodic surfaces.

## Capabilities
- **Periodicity**: 1D (nanowires/polymers), 2D (surfaces/slabs), and 3D (bulk).
- **Spectroscopy**:
  - Band structures and DOS.
  - COOP/COHP (Crystal Orbital Overlap Population) for bonding analysis.
  - EELS, STM images, Phonons.
- **Transition States**: NEB for periodic reactions.

## Key Strengths
- **Basis Set**: STOs describe electron cusps correctly and offer chemical intuition.
- **Dimensionality**: No need for artificial vacuum padding in 1D/2D systems (true 2D periodicity).
- **Analysis**: Powerful bonding analysis tools (COOP/COHP) standard in quantum chemistry but rarer in physics codes.
- **Heavy Elements**: Excellent relativistic treatment (ZORA).

## Inputs & Outputs
- **Interface**: Highly integrated with the AMS-GUI for setup and visualization.
- **Scripting**: PLAMS (Python Library for AMS) for workflow automation.
## Performance Characteristics
- **Scaling**: Formally $O(N^3)$ but utilizes linear scaling techniques for large unit cells.
- **Parallelization**: Supports hybrid MPI/OpenMP. Pure MPI is often recommended for multi-node runs; OpenMP for single node.
- **Memory**: Basis set dependent. STOs require fewer functions than GTOs for same accuracy but integral evaluation is more complex.

## Computational Cost
- **Functional Dependence**: Meta-GGAs (like SCAN) are significantly more expensive than LDAs/GGAs.
- **Relativity**: Spin-Orbit Coupling (SOC) increases cost by 4-8x compared to scalar relativistic runs.

## Limitations & Known Constraints
- **Basis Set Convergence**: High accuracy requires large basis sets (TZP, QZ4P), which can be computationally heavy.
- **Basis Linear Dependence**: Large basis sets on dense structures can lead to linear dependence issues.

## Best Practices
- **Basis Sets**: Start with a smaller basis (e.g., DZ or TZP) for geometry relaxation, then switch to larger basis for final properties.
- **Convergers**: Use `NumericalAccuracy` to fix SCF issues rather than just increasing iterations.
- **Relativity**: Use ZORA (Scalar) as default; only enable SOC if explicitly studying magnetism or heavy element splitting.

## Community and Support
- **Support**: Professional support via SCM (for license holders).
- **Tutorials**: Extensive documentation and tutorials available on the SCM website.
## Comparison with Other Codes
- **vs VASP/QE**: Plane-wave codes are often faster for simple bulk metals. BAND excels for open structures, 1D/2D systems, and systems requiring localized chemical analysis.
- **vs CRYSTAL**: Both use localized basis sets. BAND uses STOs (better cusp behavior), while CRYSTAL uses GTOs (Gaussian).

## Verification & Sources
**Primary sources**:
1.  **Official Website**: [SCM BAND](https://www.scm.com/product/band/)
2.  **Documentation**: [BAND Manual](https://www.scm.com/doc/BAND/)

**Verification status**: ✅ VERIFIED
