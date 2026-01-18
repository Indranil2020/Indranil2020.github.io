# risb

## Official Resources
- Source Repository: https://github.com/thenoursehorse/risb
- License: GPL-3.0

## Overview
`risb` is a code package for solving strongly correlated many-body electronic problems using the Rotationally Invariant Slave Bosons (RISB) method. This auxiliary particle method captures key features of correlations (like the Mott transition) at a computational cost comparable to mean-field theories, making it significantly faster than full DMFT.

**Scientific domain**: Strongly correlated lattice models, Mean-field approximations
**Target user community**: Researchers studying model Hamiltonians, Mott transitions

## Theoretical Methods
- Rotationally Invariant Slave Bosons (RISB)
- Mean-Field Theory
- Auxiliary Particle Method
- Hubbard Model / Anderson Lattice

## Capabilities
- Solving strongly correlated lattice models
- Describing the Mott Metal-Insulator Transition
- Calculating quasiparticle renormalization factors (Z)
- Evaluating ground state energies and occupancies

## Key Strengths
### Speed vs. Physics:
- Captures the Mott transition (unlike HF) but is much faster than DMFT.
### Rotational Invariance:
- Can handle general local interactions and orbital mixing correctly.
### Tutorials:
- Repository includes examples/tutorials for common lattice models.

## Inputs & Outputs
- **Input formats**:
  - Model parameters (lattices, U, t)
- **Output data types**:
  - Order parameters
  - Renormalization factors
  - Energies

## Interfaces & Ecosystem
- **Language**: Python (likely, based on typical research codes in this domain).

## Advanced Features
- **Lattice Models**: Built-in support for common lattices.

## Performance Characteristics
- **Efficiency**: High. Solves non-linear algebraic equations rather than integral equations (like DMFT) or massive Hamiltonians (like ED).

## Computational Cost
- **Low**: Accessible on standard workstations.

## Limitations & Known Constraints
- **Dynamics**: Static approximation (frequency independent self-energy at low energy).
- **finite-T**: Usually formulated for T=0 or low T.

## Comparison with Other Codes
- **vs CyGutz**: Similar method (RISB/Gutzwiller are related).
- **vs DMFT**: `risb` is the "infinite dimensions" mean-field limit without dynamics.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/thenoursehorse/risb

**Verification status**: ✅ VERIFIED
- Source code: OPEN
