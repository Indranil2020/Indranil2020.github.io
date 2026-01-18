# NCA_Standalone (CorentinB78)

## Official Resources
- Source Repository: https://github.com/CorentinB78/NCA
- License: MIT License
- Language: Python

## Overview
This is a standalone implementation of the Non-Crossing Approximation (NCA) for solving quantum impurity problems, written in Python. Unlike the TRIQS-based version, this is a self-contained code, suitable for learning or specific lightweight applications without the full TRIQS dependency stack.

**Scientific domain**: Quantum Impurity Solvers, Many-body Physics
**Target user community**: Students, Researchers needing a light, pure-Python NCA solver

## Theoretical Methods
- Non-Crossing Approximation (NCA)
- Resolvent Operator formalism
- Pseudo-particle Green's functions
- Integral equations for self-consistency

## Capabilities
- **Impurity Solver**: Solves the Single Impurity Anderson Model (SIAM).
- **Spectral Functions**: Calculates impurity spectral functions $A(\omega)$.
- **Finite Parameters**: Handles finite interaction $U$, impurity level $\epsilon_d$, and finite temperatures.
- **Hybridization**: Accepts arbitrary hybridization functions $\Delta(\omega)$.

## Key Strengths
### Simplicity:
- Pure Python implementation effectively lowers the barrier to entry.
- Minimal dependencies (NumPy/SciPy).
### Education:
- Clear code structure for understanding the integral equations of NCA.

## Inputs & Outputs
- **Input formats**:
  - Python script configuration: `U`, `ed`, `beta`, `Gamma` (hybridization width/function).
- **Output data types**:
  - Spectral densities (text files or arrays).
  - Occupancies.

## Interfaces & Ecosystem
- **Standalone**: Independent of heavy frameworks like TRIQS.
- **Integration**: Can be imported as a Python module.

## Performance Characteristics
- **Efficiency**: NCA is computationally efficient compared to QMC, but slower than simple approximations like Hubbard-I.
- **Cost**: Low execution time for typical single-impurity problems.

## Limitations & Known Constraints
- **Approximation**: NCA validity is limited (good for large N, high T, strong coupling; fails at low T Fermi liquid regime).
- **Features**: Lacks multi-orbital support and advanced features of production codes.

## Comparison with Other Codes
- **vs TRIQS-NCA**: Lighter, standalone, easier to install.
- **vs ED**: Continuous bath support, but approximate.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/CorentinB78/NCA

**Verification status**: ✅ VERIFIED
- Source code: OPEN (Python)
- Method: NCA
