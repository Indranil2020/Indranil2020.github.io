# TRIQS-cthyb

## Official Resources
- Homepage: https://triqs.github.io/cthyb/
- Documentation: https://triqs.github.io/cthyb/latest/
- Source Repository: https://github.com/TRIQS/cthyb
- License: GNU General Public License v3.0

## Overview
TRIQS/cthyb is a state-of-the-art continuous-time hybridization expansion quantum Monte Carlo impurity solver for multi-orbital Anderson impurity models. It is one of the most widely used CTQMC solvers in the DMFT community, offering efficient algorithms for solving quantum impurity problems with general multi-orbital interactions. Part of the TRIQS ecosystem, it integrates seamlessly with TRIQS/DFTTools for DFT+DMFT calculations.

**Scientific domain**: Quantum impurity problems, DMFT, strongly correlated systems  
**Target user community**: Researchers performing DMFT and DFT+DMFT calculations

## Theoretical Methods
- Continuous-time quantum Monte Carlo (CTQMC)
- Hybridization expansion (CT-HYB)
- Matrix formulation and segment picture
- Multi-orbital Anderson impurity model
- General two-body interactions (density-density and beyond)
- Good quantum numbers (particle number, spin)
- Imaginary time and Matsubara frequency formulations
- Two-particle Green's functions and vertices
- Dynamical spin and charge susceptibilities

## Capabilities (CRITICAL)
- Multi-orbital impurity problems (tested up to 5+ orbitals)
- General multi-orbital interactions (full Coulomb tensor)
- Density-density and non-density-density interactions
- Spin-orbit coupling effects
- Complex hybridization functions
- Particle-hole symmetric and asymmetric problems
- Temperature-dependent calculations
- Single-particle Green's functions and self-energies
- Two-particle correlation functions
- Improved estimators for reduced noise
- Measurement of high-frequency tails
- MPI parallelization
- Checkpoint and restart capability
- Integration with TRIQS ecosystem

**Sources**: Official TRIQS/cthyb documentation (https://github.com/TRIQS/cthyb), P. Seth et al., Comput. Phys. Commun. 200, 274 (2016), confirmed in 7/7 source lists

## Inputs & Outputs
**Input formats**:
- Python-based problem definition
- HDF5 hybridization functions
- Interaction parameters (U, J matrices)
- TRIQS Green's function objects

**Output data types**:
- Single-particle Green's functions (imaginary time and Matsubara)
- Self-energies
- Occupation numbers and double occupancies
- Two-particle Green's functions
- Monte Carlo statistics and histories
- HDF5 archives

## Interfaces & Ecosystem
- **TRIQS framework**: Native integration with TRIQS libraries
- **DFT+DMFT**: Works with TRIQS/DFTTools for ab-initio calculations
- **solid_dmft**: Used as primary impurity solver
- **Python interface**: Convenient scripting and automation
- **Analysis tools**: TRIQS-based post-processing

## Limitations & Known Constraints
- CTQMC computational cost scales with inverse temperature
- Sign problem minimal for moderate U but can appear
- Statistical errors require sufficient Monte Carlo sampling
- Multi-orbital problems memory intensive
- Two-particle quantities expensive to measure accurately
- Requires TRIQS ecosystem installation
- Learning curve for TRIQS framework


## Performance Characteristics
- **Efficiency**: State-of-the-art C++ implementation with optimized local updates.
- **Parallelization**: MPI parallelization over Monte Carlo walkers; near-linear scaling.
- **Memory**: Dense matrix operations can be memory intensive for 5+ orbitals.
- **Bottlenecks**: Matrix multiplications (BLAS level 3) and measuring two-particle quantities.

## Comparison with Other Solvers
- **vs iQIST**: iQIST offers more solver variants (CT-INT, CT-AUX) and is standalone Fortran; TRIQS/cthyb is C++/Python integrated.
- **vs w2dynamics**: Both are top-tier CT-HYB solvers; TRIQS/cthyb integrates deeply with the TRIQS library ecosystem.
- **vs ALPS/cthyb**: TRIQS/cthyb is the modern successor with better performance and active development.

## Verification & Sources
**Primary sources**:
1. Official documentation: https://triqs.github.io/cthyb/latest/
2. GitHub repository: https://github.com/TRIQS/cthyb
3. P. Seth et al., Comput. Phys. Commun. 200, 274-284 (2016) - TRIQS/cthyb paper
4. E. Gull et al., Rev. Mod. Phys. 83, 349 (2011) - CT-HYB review

**Secondary sources**:
1. TRIQS tutorials and documentation
2. Published DFT+DMFT studies using TRIQS/cthyb
3. Benchmark comparisons
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, GPL v3)
- Community support: Active (TRIQS project)
- Academic citations: >200 (main paper)
- Maintained by Flatiron Institute
- Actively developed and maintained
