# Dual Fermions (opendf / DFermion)

## Official Resources
- **Source Repository (opendf)**: https://github.com/CQMP/opendf
- **Source Repository (DFermion)**: https://github.com/huangli712/DFermion
- **Paper**: https://arxiv.org/abs/1507.00895
- **License**: GPL v3 (opendf)

## Overview
"Dual Fermions" refers to a theoretical method (Dual Fermion expansion) for going beyond Dynamical Mean Field Theory (DMFT) to include non-local spatial correlations. While not a single software package, there are specific open-source implementations, most notably **opendf** (developed by the Center for Quantum Materials Physics) and **DFermion**. These codes solve strongly correlated lattice problems (like the Hubbard model) using the dual fermion diagrammatic expansion.

**Scientific domain**: Many-body physics, Strongly Correlated Systems, DMFT extensions  
**Target user community**: Theoretical physicists

## Capabilities (CRITICAL)
- **Method**: Dual Fermion (DF) diagrammatic expansion.
- **Models**: Single-orbital Hubbard model (opendf).
- **Solvers**: Uses Continuous-Time Quantum Monte Carlo (CT-QMC) as an impurity solver (often external or integrated).
- **Correlations**: Captures non-local spatial correlations missed by standard DMFT.
- **Dimensions**: Applicable to 2D square, 3D cubic lattices, etc.

## Inputs & Outputs
- **Input formats**: Configuration files (parameters for interaction U, temperature T, chemical potential).
- **Output data types**: Self-energies, Green's functions, susceptibility.

## Performance Characteristics
- Computationally intensive due to Monte Carlo sampling and diagrammatic summation.
- Scales with the number of diagrams/orders considered.

## Application Areas
- Phase transitions in Hubbard models (e.g., antiferromagnetism).
- Pseudogap physics in cuprates.
- Non-local correlation effects.

## Verification & Sources
- **opendf**: Verified existence on GitHub (CQMP/opendf) and published paper (Phys. Procedia).
- **DFermion**: Verified existence on GitHub.

**Confidence**: VERIFIED (as specific implementations of the method)

**Verification status**: ✅ VERIFIED
- Website: GITHUB REPO
- Documentation: PAPER / README
- Source: OPEN
- Development: RESEARCH CODE (Niche)
- Applications: Dual Fermion Method
