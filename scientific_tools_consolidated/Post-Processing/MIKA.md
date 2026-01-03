# MIKA

## Official Resources
- Homepage: https://wiki.fysik.dtu.dk/mika/
- Documentation: https://wiki.fysik.dtu.dk/mika/
- Source Repository: Distributed via website (GPL)
- License: GNU General Public License

## Overview
MIKA (Multigrid Instead of K-spAce) is a collection of Matlab/Octave functions and C++ codes for electronic structure calculations using real-space grid methods (finite difference and multigrid). It includes a DFT solver (RMG) and a time-dependent DFT solver. It is designed for solving the Schrödinger and Poisson equations on real-space grids, particularly useful for transport and large systems without periodic boundary conditions.

**Scientific domain**: Real-space DFT, multigrid methods, electronic structure  
**Target user community**: Developers of real-space methods, researchers in transport

## Theoretical Methods
- Real-space finite difference discretization
- Multigrid acceleration for Poisson/Schrödinger solvers
- Density Functional Theory (DFT)
- Time-Dependent DFT (TDDFT)
- Generalized Poisson equation

## Capabilities (CRITICAL)
- **RMG**: Real-space Multigrid DFT code (included in MIKA suite)
- **FD**: Finite difference derivatives
- **Poisson Solver**: Fast multigrid Poisson solver
- **Transport**: Wavefunction matching for transport calculations (Do calculations for scattering states)
- **Jellium**: Jellium model calculations

**Sources**: MIKA website, Comp. Phys. Comm. 128, 1 (2000)

## Inputs & Outputs
- **Input formats**: Matlab/Octave scripts, structure files
- **Output data types**: Grid data (potential, density), wavefunctions

## Interfaces & Ecosystem
- **Matlab/Octave**: Core environment
- **C++**: Performance-critical parts
- **GPAW**: Similar real-space philosophy (developed by related groups at DTU)

## Workflow and Usage
1. Define grid and potential in Matlab.
2. Call MIKA solvers (e.g., `mg_solve`).
3. Analyze wavefunctions on the grid.

## Performance Characteristics
- O(N) scaling for Poisson solver
- Efficient for non-periodic systems (clusters, wires)
- Limited by grid memory usage

## Application Areas
- Nanowire transport
- Quantum dots
- Method development for real-space DFT
- Jellium clusters

## Community and Support
- Developed at CSC (Finland) and DTU (Denmark)
- Legacy code (precursor to modern real-space codes like GPAW/Octopus)
- Open-source

## Verification & Sources
**Primary sources**:
1. Homepage: https://wiki.fysik.dtu.dk/mika/
2. Publication: T. Torsti et al., Comp. Phys. Comm. 162, 167 (2004)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (Wiki)
- Documentation: AVAILABLE
- Source: OPEN (GPL)
- Development: STABLE (Legacy/Maintenance)
- Applications: Real-space DFT, multigrid
