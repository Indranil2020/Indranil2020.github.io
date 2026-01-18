# merapp

## Official Resources
- **Repository**: https://github.com/g1257/merapp
- **License**: GNU General Public License v3.0
- **Developer**: Gonzalo Alvarez (g1257)

## Overview
merapp (MERA++) is a robust C++ implementation of the Multi-scale Entanglement Renormalization Ansatz (MERA) algorithm. It is part of a suite of tools (including DMRG++) designed for the simulation of strongly correlated quantum systems. MERA++ is specifically engineered to handle scale-invariant systems and quantum critical points by optimizing the MERA tensor network, which adds an extra dimension of "scale" to efficiently capture critical entanglement.

**Scientific domain**: Quantum Critical Phenomena, Conformal Field Theory (CFT), Condensed Matter.
**Target user community**: Physicists studying quantum phase transitions and critical points.

## Theoretical Methods
- **MERA**: Ternary, Binary, and Modified Binary MERA ansatzes.
- **Entanglement Renormalization**: Using isometries and disentanglers to coarse-grain the lattice.
- **Variational Optimization**: Minimizing energy expectation value.
- **Scale Invariance**: Determining scaling dimensions and central charge.

## Capabilities (CRITICAL)
- **Generic Engine**: Supports arbitrary dimensions, arities, and geometries.
- **Hamiltonians**: Built-in support for Heisenberg, Hubbard, t-J, and other models.
- **Symmetries**: U(1), Z2, and SU(2) symmetry support via `PsimagLite` backend.
- **Observables**: Computation of energy, correlation functions, and scaling operators.
- **Parallelism**: Shared memory parallelization (pthreads/OpenMP).

## Key Strengths
### Specialization
- One of the few open-source, production-ready codes specifically for MERA.
- Handles the complex geometry of MERA networks efficiently.

### Performance
- Written in optimized C++.
- Exploits symmetries to drastically reduce tensor sizes.

## Inputs & Outputs
- **Input**: Input file defining model parameters, MERA type, and symmetry sector.
- **Output**: Ground state energy, Optimized tensors, Conformal data (scaling dims).

## Interfaces & Ecosystem
- **Dependencies**: `PsimagLite` (linear algebra/utils library from same author).
- **Interface**: Command-line driven.
- **Ecosystem**: Compatible with `DMRG++` input formats and utilities.

## Advanced Features
- **Conformal Data**: Extraction of central charge and primary field scaling dimensions.
- **Restart**: Continuation of optimization from saved states.

## Performance Characteristics
- **Cost**: MERA contraction is $O(\chi^7)$ or similar, making it expensive but polynomial.
- **Optimization**: Uses efficient local update strategies.

## Computational Cost
- **High**: Significantly more expensive than DMRG for the same bond dimension, but captures criticality better.

## Comparison with Other Codes
- **vs DMRG**: DMRG (MPS) fails to capture critical entanglement (logarithmic growth) efficiently; MERA succeeds but at higher constant cost.
- **vs TeNPy**: TeNPy focuses on MPS; MERA++ focuses on MERA.

## Application Areas
- **Quantum Critical Points**: Ising, Potts, Heisenberg transitions.
- **Topological Phases**: Study of topological order using MERA.

## Best Practices
- **Symmetries**: Use symmetries whenever possible to enable larger $\chi$.
- **Convergence**: Monitor energy and scaling dimensions for stability.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/g1257/merapp
2. Associated papers by G. Alvarez on MERA/DMRG.

**Confidence**: VERIFIED - Open source distribution.
**Verification status**: ✅ VERIFIED
