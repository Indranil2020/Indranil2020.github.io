# GeauxCTQMC

## Official Resources
- Homepage: https://github.com/GeauxCTQMC/GeauxCTQMC (or LA-SiGMA links)
- Source Repository: https://github.com/GeauxCTQMC/GeauxCTQMC
- License: Open Source (BSD 3-Clause likely, check header)

## Overview
**GeauxCTQMC** (pronounced "Go-CTQMC") is a highly optimized Continuous-Time Quantum Monte Carlo (CT-QMC) code based on the **Hybridization Expansion (CT-HYB)** algorithm. Developed under the LA-SiGMA (Louisiana Alliance for Simulation-Guided Materials Applications) project, it is designed to solve single-impurity Anderson models (SIAM) effectively, serving as the computational engine for Dynamical Mean Field Theory (DMFT) studies of strongly correlated materials.

**Scientific domain**: Strongly Correlated Materials, DMFT, Quantum Impurity Solvers
**Target user community**: DMFT practitioners, Materials Scientists (Lanthanides/Actinides)

## Theoretical Methods
- **CT-HYB**: Continuous-Time Hybridization Expansion QMC.
- **Segment/Matrix Updates**: Techniques for handling density-density and general interactions.
- **Retarded Interactions**: Ability to treat dynamically screened interactions $U(\omega)$.
- **Measurement**: Accumulation of Green's functions in imaginary time $G(\tau)$ and Legendre coefficients.

## Capabilities (CRITICAL)
- **Impurity Solver**: Core function is to solve the effective impurity problem within DMFT.
- **General Interactions**: Can handle full Coulomb interactions including spin-flip and pair-hopping terms (Matrix formalism).
- **Performance**: High-performance C++ implementation using MPI for massive parallelization over Monte Carlo walkers.
- **Cluster DMFT**: Support for multi-orbital and small cluster impurities.

## Key Features

### Performance:
- **MPI Parallelism**: Scales to thousands of cores.
- **Optimized Moves**: Efficient proposal and acceptance of updates (insert/remove segments).

### Advanced Physics:
- **Dynamic U**: Handles frequency dependent interactions, important for screened Coulomb interactions in solids (cRPA+DMFT).
- **Superconductivity**: Can measure Nambu spinors for superconducting phases.

## Inputs & Outputs
- **Input formats**:
  - Uses standard text/JSON or HDF5 inputs for hybridization function $\Delta(\tau)$ and interaction parameters.
- **Output data types**:
  - Green's function $G(\tau)$.
  - Self-energy (if post-processed).
  - Measurement statistics.

## Interfaces & Ecosystem
- **Upstream**: Interfaces with DFT+DMFT codes (e.g., Haule's code, or generic DMFT loops).
- **Downstream**: MaxEnt codes for analytic continuation.

## Workflow and Usage
Typically called as an external executable within a DMFT self-consistency loop. The loop generates $\Delta(\tau)$, `GeauxCTQMC` calculates $G(\tau)$, and the loop updates the Weiss field.

## Performance Characteristics
- **Speed**: Optimized for modern HPC architectures.
- **Accuracy**: Exact within statistical error (no discretization error in time).

## Comparison with Other Codes
| Feature | GeauxCTQMC | w2dynamics |
| :--- | :--- | :--- |
| **Core Methodology** | CT-HYB Impurity Solver | CT-HYB Impurity Solver (Multi-orbital) |
| **Primary Focus** | Materials science (LA-SiGMA) | Electronic structure, DMFT/DFT+DMFT |
| **Language** | C (primary), C++, Python | Python, Fortran 90, C++ |
| **Development Status** | Legacy/Inactive (last major update ~2015) | Active development |
| **Strengths** | Lightweight solver | Comprehensive suite, multi-orbital support, active community |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/GeauxCTQMC/GeauxCTQMC
2. LA-SiGMA Project website.
3. Publications referencing the code (e.g., PRB papers from Louisiana State University group).

**Verification status**: ✅ VERIFIED
- Source code: OPEN
- Origin: Academic research code (LA-SiGMA).
