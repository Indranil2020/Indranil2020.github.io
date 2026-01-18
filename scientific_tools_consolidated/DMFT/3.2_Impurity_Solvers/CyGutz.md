# CyGutz

## Official Resources
- Source Repository: https://github.com/yaoyongxin/CyGutz
- License: Custom Open Source (Ames Lab / Rutgers / DOE - BSD-like)

## Overview
CyGutz is a Gutzwiller solver implemented in Cython/Python. It is designed to solve generic tight-binding models with local interactions using the Gutzwiller-Rotationally Invariant Slave-Boson (RISB) method. It optimizes the single Slater determinant and local many-body degrees of freedom simultaneously within the Gutzwiller approximation.

**Scientific domain**: Strongly correlated electrons, Gutzwiller approximation, Slave-Boson methods
**Target user community**: Researchers studying Hubbard models, tight-binding systems with correlations

## Theoretical Methods
- Gutzwiller Approximation
- Rotationally Invariant Slave-Boson (RISB) theory
- Variational Monte Carlo (connection to)
- Tight-Binding models

## Capabilities
- Solving generic tight-binding models with local interactions
- Optimizing Gutzwiller variational parameters
- Handling multi-orbital systems
- Calculation of renormalized band structures

## Key Strengths
### Efficiency:
- Mean-field cost, much cheaper than DMFT or QMC.
### Generic Models:
- Can handle general tight-binding Hamiltonians.
### Implementation:
- Cython helps in performance while maintaining Python usability.

## Inputs & Outputs
- **Input formats**:
  - Tight-binding parameters (hopping)
  - Interaction parameters (Coulomb U, J)
- **Output data types**:
  - Renormalized hoppings
  - Quasiparticle weights (Z)
  - Ground state energy
  - Orbital occupations

## Interfaces & Ecosystem
- **Python**: Designed to be used as a Python library/module.

## Advanced Features
- **Rotationally Invariant**: Handles general interaction tensors effectively.
- **Simultaneous Optimization**: Optimizes both the uncorrelated wavefunction and the projector.

## Performance Characteristics
- **Speed**: Very fast compared to dynamic solvers.
- **Scaling**: Scales reasonably with system size, limited mainly by the local Hilbert space size for the slave bosons.

## Computational Cost
- **Low**: Mean-field level cost.

## Limitations & Known Constraints
- **Approximation**: It is a static mean-field theory (Gutzwiller), missing dynamic correlations (frequency dependence).
- **Accuracy**: Good for ground state properties (Fermi surface, mass enhancement) but fails for satellites/incoherent features.

## Comparison with Other Codes
- **vs DMFT**: Cheaper, static, no spectral functions (only quasiparticles).
- **vs DFT+U**: More flexible treatment of local correlations (screening, mass enhancement).

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/yaoyongxin/CyGutz

**Verification status**: ✅ VERIFIED
- Source code: OPEN
