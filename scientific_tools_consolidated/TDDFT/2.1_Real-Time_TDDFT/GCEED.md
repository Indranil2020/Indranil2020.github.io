# GCEED

## Official Resources
- **SourceForge**: https://sourceforge.net/projects/gceed/
- **Website**: https://gceed.nict.go.jp/
- **License**: GNU General Public License (GPL)

## Overview
GCEED (Grid-based Coupled Electron and Electromagnetic field Dynamics) is an open-source software package designed for massively parallel first-principles calculations of electron dynamics in real time and real space. It is specifically built to simulate the coupled dynamics of electrons and electromagnetic fields, making it suitable for studying light-matter interactions in nanostructures and solids.

**Scientific domain**: Nanostructures, solids, intense laser-matter interactions, plasmonics
**Target user community**: Researchers simulation ultrafast electron dynamics and light-matter coupling

## Theoretical Methods
- **Time-Dependent Density Functional Theory (TDDFT)**: Real-time propagation
- **Maxwell-TDDFT**: Self-consistent solution of Maxwell's equations and Time-Dependent Kohn-Sham (TDKS) equations
- **Real-Space Grid**: Finite-difference representation for wavefunctions and fields
- **Pseudopotentials**: Norm-conserving pseudopotentials

## Capabilities
- **Couple Dynamics**: Simultaneous propagation of electron density and electromagnetic fields
- **Intense Fields**: Simulation of non-linear optical phenomena and strong field physics
- **Massive Parallelism**: Efficient scaling on supercomputers
- **Flexible Boundary Conditions**: Absorbing boundary conditions for open systems

## Inputs & Outputs
- **Input formats**:
  - Text-based input configuration files
  - Pseudopotential files (standard formats)
  - Coordinate files for atomic structures
- **Output data types**:
  - Time-dependent densities and fields
  - Induced currents and dipoles
  - High-harmonic generation (HHG) spectra
  - Near-field distributions

## Performance Characteristics
- **Parallelization**: Hybrid MPI/OpenMP implementation
- **Scaling**: Designed for large-scale calculations on supercomputers
- **Grid Efficiency**: Real-space grid avoids basis set linear dependence issues

## Computational Cost
- **Scaling**: Linear scaling with respect to system size for the grid operations, but dominated by propagation steps.
- **Memory**: Efficient distribution of grid points across processors.

## Limitations & Known Constraints
- **Basis Set**: Real-space grids require careful convergence testing of grid spacing.
- **Documentation**: Primary documentation is via the manual and SourceForge page, may be less extensive than major commercial codes.

## Comparison with Other Codes
- **vs Octopus**: Both are real-space real-time TDDFT codes; GCEED emphasizes coupled Maxwell-TDDFT dynamics for intense fields.
- **vs SALMON**: Similar capabilities for light-matter interaction; GCEED has specific optimizations for coupled field dynamics.

## Best Practices
- **Grid Spacing**: Ensure grid is fine enough to describe high-energy states if needed.
- **Time Step**: Stability of the propagation depends on the Courant condition for the Maxwell solver.

## Citations
- **Primary**: See SourceForge page for specific development team citations (NICT/IMS Japan).
