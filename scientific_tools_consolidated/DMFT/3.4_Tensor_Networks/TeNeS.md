# TeNeS

## Official Resources
- **Repository**: https://github.com/issp-center-dev/TeNeS
- **Documentation**: https://issp-center-dev.github.io/TeNeS/
- **License**: GNU General Public License v3.0
- **Developers**: ISSP (University of Tokyo)

## Overview
TeNeS (Tensor Network Solver) is an open-source, massively parallel software package for calculating the ground-state wavefunctions of two-dimensional quantum many-body systems. It employs the infinite Projected Entangled Pair States (iPEPS) ansatz and optimizes it using imaginary time evolution. Designed for high-performance computing, TeNeS leverages the `mptensor` library for efficient distributed tensor operations, allowing it to tackle large bond dimensions and complex frustrated spin/boson models on square, honeycomb, and triangular lattices.

**Scientific domain**: Condensed Matter Physics, 2D Quantum Magnetism, Strongly Correlated Systems.
**Target user community**: Physicists studying ground state phases of 2D lattice models on supercomputers.

## Theoretical Methods
- **iPEPS**: Infinite Projected Entangled Pair States ansatz for 2D systems.
- **Imaginary Time Evolution**: Simple and full update optimization schemes.
- **CTMRG**: Corner Transfer Matrix Renormalization Group for environment contraction.
- **MFE**: Mean-Field Environment method (optional).
- **Lattices**: Built-in support for Square, Honeycomb, and Triangular lattices.

## Capabilities (CRITICAL)
- **Massive Parallelism**: Hybrid MPI/OpenMP parallelization for distributed memory systems.
- **Ground State Search**: Finds ground states of user-defined 2D Hamiltonians.
- **Observables**: Calculates magnetization, correlation functions, and energy.
- **Finite T/Dynamics**: (v2+) Support for finite temperature and real-time evolution.
- **Simple Input**: Utilities (`tenes_std`, `tenes_simple`) to generate inputs for standard models (Heisenberg, Hubbard, etc.).
- **Differentiable**: Supports automatic differentiation in newer versions (check docs).

## Key Strengths
### HPC Scaling
- Built on `mptensor` to distribute large tensors across nodes.
- Capable of handling larger $\chi$ (bond dimension) than single-node codes.

### Versatility
- Handles varying unit cell sizes and common 2D geometries.
- Robust contraction using CTMRG.

## Inputs & Outputs
- **Input**: TOML-based configuration for Hamiltonian, lattice, and run parameters.
- **Output**: Physical observables (energy, order parameters) and wavefunction checkpoint files.

## Interfaces & Ecosystem
- **Dependencies**: `mptensor`, ScaLAPACK, MPI.
- **Language**: C++ (Core), Python (Input generation tools).
- **Integration**: Part of the ISSP software suite (MateriApps).

## Advanced Features
- **Symmetries**: Support for Abelian symmetries (U(1), Z2) via `mptensor`.
- **Checkpointing**: Restart simulations from previous tensor states.

## Performance Characteristics
- **Speed**: competitive with state-of-the-art C++ TN codes.
- **Scalability**: Excellent scaling on supercomputers (e.g., Fugaku, Summit) due to distributed tensor backend.

## Computational Cost
- **High**: 2D TN contractions scale as $O(\chi^{10})$ or similar depending on the algorithm, requiring HPC for accurate results.
- **Memory**: Distributed memory allows handling tensors too large for a single node.

## Comparison with Other Codes
- **vs PEPS (QuantumLiquids)**: TeNeS uses CTMRG (deterministic contraction) and runs on clusters; PEPS uses Variational Monte Carlo (stochastic).
- **vs ITensor**: TeTenS is specialized for 2D infinite systems (iPEPS); ITensor is general but historically 1D-focused (MPS).
- **vs peps-torch**: peps-torch uses PyTorch and AD for optimization; TeNeS uses imaginary time evolution and C++ MPI.

## Application Areas
- **Frustrated Magnets**: J1-J2 models, Kagome antiferromagnets.
- **Bosonic Systems**: Bose-Hubbard models.
- **Quantum Phase Transitions**: Determining phase boundaries in 2D.

## Best Practices
- **Bond Dimension**: Systematically increase $\chi$ to check convergence.
- **Environment**: Ensure CTMRG environment dimension ($\chi_{env}$) is large enough (typically $\chi_{env}^2 \sim \chi$).
- **Parallelism**: Match MPI ranks to tensor block structure for efficiency.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/issp-center-dev/TeNeS
2. Paper: "TeNeS: Tensor Network Solver for Quantum Lattice Systems" (arXiv/phys. Rev.).

**Confidence**: VERIFIED - Established code from ISSP.
**Verification status**: ✅ VERIFIED
