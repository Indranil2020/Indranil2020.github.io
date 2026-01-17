# MuST (Multiple Scattering Theory Suite)

## Official Resources
- Homepage: https://github.com/mstsuite/MuST
- Documentation: https://must.readthedocs.io/
- Source Repository: https://github.com/mstsuite/MuST
- License: BSD 3-Clause

## Overview
MuST is an open-source ab initio electronic structure suite based on Multiple Scattering Theory (MST). It integrates the Korringa-Kohn-Rostoker (KKR) Green function method and the Locally Self-consistent Multiple Scattering (LSMS) method. It is unique in its ability to handle disordered materials (via CPA) and scale to petascale/exascale computing systems for massive all-electron calculations.

**Scientific domain**: Disordered materials, alloys, magnetic systems, quantum phase transitions
**Target user community**: Researchers in metallurgy, disordered systems, and high-performance computing

## Theoretical Methods
- Density Functional Theory (DFT)
- Multiple Scattering Theory (MST) / KKR
- Coherent Potential Approximation (CPA) for random alloys
- Locally Self-consistent Multiple Scattering (LSMS)
- Kubo-Greenwood formula for conductivity
- Landau-Lifshitz-Gilbert (LLG) dynamics for spins
- Full-potential and Muffin-tin approximations
- Relativistic effects (Scalar and Fully Relativistic)

## Capabilities
- Electronic structure of ordered and disordered solids
- Linear-scaling (O(N)) calculations for tens of thousands of atoms
- Electrical conductivity in random alloys
- Spin dynamics and thermodynamics
- First-principles calculation of critical temperatures
- Defect states and impurities

## Key Strengths
### Disordered Materials
- KKR-CPA method efficiently averages over disorder without supercells
- Direct calculation of alloy properties

### Massive Scaling (LSMS)
- O(N) scaling allows simulation of extremely large systems
- Designed for top-tier supercomputers (Exascale ready)
- GPU acceleration

## Inputs & Outputs
- **Input**: Fortran-namelist style inputs, atom positions, potential files
- **Output**: DOS, Band structure (KKR), conductivity, magnetic moments

## Interfaces & Ecosystem
- **Input/Output**:
  - Fortran namelist input format.
  - HDF5 support for large data handling.
  - Integration with multiple scattering theory analysis tools.
- **Libraries**:
  - Uses LAPACK/BLAS, MPI, and OpenMP.

## Advanced Features
- **KKR-CPA**:
  - Coherent Potential Approximation for rigorous treatment of random alloys.
  - Handles chemical disorder without large supercells.
- **Kubo-Greenwood Transport**:
  - First-principles calculation of electrical conductivity.
  - Residual resistivity in alloys.
- **Landau-Lifshitz-Gilbert (LLG)**:
  - Spin dynamics simulations for magnetic systems.
  - Thermodynamic properties of magnets.

## Community and Support
- **Documentation**: Comprehensive guides at https://must.readthedocs.io/
- **Development**: Managed by Oak Ridge National Lab (ORNL) and collaborators.
- **Support**: GitHub issues and documentation tutorials.

## Computational Cost
- **LSMS**: Linear scaling O(N), highly efficient for very large supercells.
- **KKR**: Efficient for periodic unit cells and alloys using Green's functions.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/mstsuite/MuST
2. Documentation: https://must.readthedocs.io/
3. "MuST: An open source package for ab initio electronic structure calculations in diverse computing environments"

**Confidence**: VERIFIED
**Status**: Open Source, Active (ORNL/NSF supported)
