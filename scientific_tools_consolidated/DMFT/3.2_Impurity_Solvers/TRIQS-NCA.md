# TRIQS-NCA

## Official Resources
- Source Repository: https://github.com/amoutenet/NCA
- License: Open Source (Check repository)

## Overview
TRIQS-NCA is an implementation of the Non-Crossing Approximation (NCA) impurity solver, developed to work within the TRIQS (Toolbox for Research on Interacting Quantum Systems) ecosystem. It provides a diagrammatic solver for the Anderson Impurity Model features, effective for problems where the hybridization is the perturbation.

**Scientific domain**: Quantum Impurity Solvers, Diagrammatic Monte Carlo/Approximations
**Target user community**: TRIQS users, researchers needing diagrammatic solvers

## Theoretical Methods
- Non-Crossing Approximation (NCA)
- Diagrammatic Expansion (lowest order in 1/N or hybridization)
- Green's Functions

## Capabilities
- Solving the Anderson Impurity Model
- Working within the TRIQS framework
- Calculating self-energies and Green's functions
- Good for high degeneracy / large Coulomb repulsion regimes (Kondo/mixed valence)

## Key Strengths
### TRIQS Integration:
- Directly compatible with the TRIQS library, allowing use in larger DMFT loops driven by TRIQS.
### Physics Regime:
- Efficient for capturing Kondo features and atomic limit physics where some other solvers struggle or are expensive.
### Spectral Resolution:
- Directly works on the real axis (or close to it) in some formulations, or imaginary axis. (NCA is often formulated on Real axis).

## Inputs & Outputs
- **Input formats**:
  - TRIQS Green's function objects
  - Hamiltonian/Interaction parameters
- **Output data types**:
  - Solved Green's function (G)
  - Self-energy (Sigma)

## Interfaces & Ecosystem
- **Dependency**: Requires TRIQS library.
- **Language**: C++ / Python (TRIQS standard).

## Advanced Features
- **Collaboration**: Developed with experts in the field (Georges, Parcollet et al.).

## Performance Characteristics
- **Speed**: Generally faster than CT-QMC, especially at low temperatures vs high frequencies.
- **Accuracy**: Approximate. Fails at very low temperatures (below Kondo temperature artifacts) in some formulations.

## Computational Cost
- **Moderate**: Cheaper than exact QMC.

## Limitations & Known Constraints
- **Approximation**: NCA is not exact. It misses vertex corrections (crossing diagrams).
- **Artifacts**: Can show non-Fermi liquid artifacts at T -> 0 in single channel models.

## Comparison with Other Codes
- **vs CT-HYB**: NCA is approximate, CT-HYB is exact. NCA is faster and has no sign problem.
- **vs OCA**: One Crossing Approximation is the next order improvement over NCA.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/amoutenet/NCA

**Verification status**: ✅ VERIFIED
- Source code: OPEN
