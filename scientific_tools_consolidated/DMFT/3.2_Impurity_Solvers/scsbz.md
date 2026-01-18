# scsbz (Slave-Boson Superconductivity)

## Official Resources
- Source Repository: https://github.com/tflovorn/scsbz
- License: Open Source (Check repository)

## Overview
`scsbz` is a code designed to solve the mean-field self-consistent equations for slave-boson superconductivity. It is based on the Kotliar-Liu slave-boson formalism (PRB 38, 7, 1988), typically applied to the study of high-temperature superconductivity in cuprate materials.

**Scientific domain**: Unconventional Superconductivity, Slave-Boson Theory
**Target user community**: Researchers in strongly correlated superconductivity

## Theoretical Methods
- Slave-Boson Mean-Field Theory
- Kotliar-Liu Formalism
- Superconductivity (d-wave, etc.)
- Hubbard / t-J Models

## Capabilities
- Solving self-consistent mean-field equations
- Describing the superconducting order parameter in correlated systems
- Calculating critical temperatures and phase diagrams (in mean-field)

## Key Strengths
### Specificity:
- Tailored for the slave-boson description of superconductivity.
### Implementation:
- Provides a reference implementation for these specific non-linear equations.

## Inputs & Outputs
- **Input formats**:
  - Model parameters (t, J, doping)
- **Output data types**:
  - Order parameters (gap, boson fields)
  - Free energies

## Interfaces & Ecosystem
- **Language**: C++ / Python (Inferred).

## Advanced Features
- **Symmetry breaking**: Can handle superconducting order parameters.

## Performance Characteristics
- **Efficiency**: Mean-field cost, efficient.

## Computational Cost
- **Low**: Standard workstation efficiency.

## Limitations & Known Constraints
- **Mean-Field**: Neglects fluctuations (gauge fields).
- **Model specific**: Tailored to cuprate-like models.

## Comparison with Other Codes
- **vs DMFT**: Mean-field slave boson is a static approximation compared to dynamic DMFT.
- **vs Bogoliubov-de Gennes**: Includes correlation effects via slave bosons unlike standard BdG.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/tflovorn/scsbz

**Verification status**: ✅ VERIFIED
- Source code: OPEN
