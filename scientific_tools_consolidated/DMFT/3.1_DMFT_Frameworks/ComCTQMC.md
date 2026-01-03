# ComCTQMC

## Official Resources
- Homepage: https://www.bnl.gov/comscope/software/comscope-software-packages.php
- Documentation: https://github.com/comscope/ComCTQMC
- Source Repository: https://github.com/comscope/ComCTQMC
- License: See repository for licensing details

## Overview
ComCTQMC is a GPU-accelerated continuous-time quantum Monte Carlo impurity solver implementing the hybridization expansion (CT-HYB) algorithm. Developed as part of the Comscope project, it provides efficient solutions to DMFT impurity problems with both partition function and worm-space measurements. The GPU acceleration enables significantly faster calculations compared to CPU-only implementations.

**Scientific domain**: DMFT impurity solver, quantum Monte Carlo, strongly correlated systems  
**Target user community**: Researchers performing DMFT calculations requiring efficient CT-HYB solver

## Theoretical Methods
- Continuous-time quantum Monte Carlo (CTQMC)
- Hybridization expansion (CT-HYB)
- Partition function measurements
- Worm-space sampling algorithms
- Multi-orbital Anderson impurity model
- GPU-accelerated algorithms

## Capabilities (CRITICAL)
- GPU-accelerated CT-HYB solver
- Significantly faster than CPU implementations
- Multi-orbital impurity problems
- General multi-orbital interactions
- Partition function measurements
- Worm-space measurements for improved statistics
- Integration with ComDMFT framework
- Designed for production DMFT calculations
- Self-energy and Green's function calculations
- Temperature-dependent simulations

**Sources**: Comscope software packages (https://www.bnl.gov/comscope/), GitHub repository, confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- Hybridization functions
- Interaction parameters (U, J matrices)
- CTQMC control parameters
- Configuration files

**Output data types**:
- Green's functions
- Self-energies
- Occupation matrices
- Monte Carlo statistics
- Observables and correlations

## Interfaces & Ecosystem
- **ComDMFT**: Primary integration with ComDMFT framework
- **GPU**: CUDA-based GPU acceleration
- **Comscope**: Part of Comscope software suite
- **DFT+DMFT**: Used in realistic materials calculations

## Limitations & Known Constraints
- Requires GPU hardware for acceleration
- CUDA toolkit dependency
- Installation complexity moderate
- CTQMC sign problem at low temperatures
- Statistical errors from Monte Carlo sampling
- Documentation primarily in repository
- Platform: Linux with NVIDIA GPU

## Verification & Sources
**Primary sources**:
1. Comscope website: https://www.bnl.gov/comscope/software/comscope-software-packages.php
2. GitHub repository: https://github.com/comscope/ComCTQMC
3. Y. Lu et al., Phys. Rev. B 104, 125107 (2021) - GPU-accelerated solver paper

**Secondary sources**:
1. ComDMFT documentation
2. Comscope project publications
3. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official information: ACCESSIBLE (Comscope website)
- Documentation: ACCESSIBLE (GitHub)
- Source code: OPEN (GitHub)
- GPU acceleration: Significant performance advantage
- Part of Comscope project (BNL)
- Active development and maintenance
