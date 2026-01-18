# CTAUX

## Official Resources
- Homepage: https://github.com/danielguterding/ctaux
- Source Repository: https://github.com/danielguterding/ctaux
- License: GPL-3.0

## Overview
**CTAUX** is a Continuous-Time Auxiliary Field Quantum Monte Carlo solver for quantum impurity models. Developed largely by Daniel Guterding, it implements the **CT-AUX** algorithm, which expands the partition function in powers of the interaction strength $U$. This method is complementary to CT-HYB and is particularly efficient for models with large coordination numbers or specific interaction forms where the weak-coupling expansion converges rapidly.

**Scientific domain**: Computational Many-Body Physics, DMFT
**Target user community**: Researchers studying Hubbard models, Cluster DMFT

## Theoretical Methods
- **CT-AUX**: Continuous-Time Auxiliary Field expansion.
- **Legendre Basis**: Representations of Green's functions using Legendre polynomials for compact storage and noise reduction.
- **Nambu Formalism**: Implementation often supports Nambu spinor formulation for superconductivity.
- **Sub-Matrix Updates**: Fast update algorithms for determinants.

## Capabilities (CRITICAL)
- **Cluster DMFT**: Well-suited for DCA (Dynamical Cluster Approximation) or CDMFT (Cellular DMFT) due to scaling properties.
- **Multi-orbital**: Handles multi-orbital impurities effectively.
- **Observables**: One- and two-particle Green's functions, spin/charge susceptibilities.

## Key Features

### Efficiency:
- **Legendre Polynomials**: Significant memory and I/O reduction by storing G in Legendre basis instead of dense time grids.
- **MPI Parallelism**: Distributes walkers across cores.

### Usability:
- **C++**: Modern C++ implementation.
- **Integration**: Can be linked with GSL (GNU Scientific Library).

## Inputs & Outputs
- **Input formats**:
  - Parameter files specifying $U$, $\beta$, chemical potential, and bath Green's function $G_0$.
- **Output data types**:
  - Interacting Green's function $G_{imp}$ in Legendre or time representation.

## Interfaces & Ecosystem
- **Dependencies**: GSL, MPI, BLAS/LAPACK.
- **Context**: Often used in conjunction with simplified DMFT loops or specific research projects on superconductivity.

## Workflow and Usage
Standard DMFT solver usage:
1. Read $G_0(\tau)$.
2. Run Monte Carlo updates (Insert/Remove auxiliary spins).
3. Measure $G(\tau)$.
4. Output results.

## Performance Characteristics
- **Regime**: Superior to CT-HYB at weak-to-intermediate interactions.
- **Scaling**: $O(k^3)$ where $k$ is expansion order (related to $U \beta$).

## Comparison with Other Codes
| Feature | CTAUX (Guterding) | CT-HYB (TRIQS) | w2dynamics |
| :--- | :--- | :--- | :--- |
| **Algorithm** | CT-AUX | CT-HYB | CT-HYB / CT-INT |
| **Coupling Regime** | Weak / Intermediate | Strong | All |
| **Implementation** | Standalone C++ | C++ / Python Lib | C++ / Python |
| **Specialization** | Clusters / Legendre Basis | General Impurity | Multi-orbital / Clusters |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/danielguterding/ctaux
2. D. Guterding's academic publications.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GPL-3.0)
- Method: Standard CT-AUX algorithm.
