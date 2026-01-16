# cmpy (Condensed Matter Python)

## Official Resources
- **Homepage**: https://github.com/dylanljones/cmpy
- **Source Repository**: https://github.com/dylanljones/cmpy
- **License**: MIT License
- **PyPI**: Not listed (install via git)

## Overview
**cmpy** is a collection of Python tools designed for computational condensed matter physics. It provides a flexible framework for constructing and solving many-body Hamiltonians, specifically focusing on lattice models like the Hubbard model and Anderson impurity model. The library includes utilities for exact diagonalization, Green's function calculations, and tight-binding models.

**Scientific domain**: Condensed Matter Physics, Lattice Models, Many-Body Physics
**Target user community**: Researchers and students working on lattice models, exact diagonalization, and Green's functions.

## Theoretical Methods
- **Exact Diagonalization (ED)**: For solving small many-body clusters.
- **Green's Functions**: Tools for calculating and manipulating Green's functions.
- **Tight-Binding Models**: Construction of lattice Hamiltonians.
- **Hubbard Model**: Specific implementations for Hubbard interactions.
- **Anderson Impurity Model**: Tools for impurity problems.
- **Linear Operators**: Efficient handling of sparse Hamiltonian matrices.

## Capabilities
- **Basis Construction**: Automatic generation of many-body basis states for fermions.
- **Sector Decomposition**: Handling of symmetry sectors (particle number, spin).
- **Hamiltonian Construction**:
  - Hopping terms
  - Interaction terms (Hubbard U)
  - Sparse matrix representation (via SciPy)
- **State Analysis**:
  - Occupation numbers
  - Spin states
  - Basis state labeling
- **Interoperability**: Built on top of NumPy and SciPy.

## Key Strengths
- **Python-Native**: Fully written in Python for ease of use and modification.
- **Educational**: Clear structure suitable for learning many-body physics concepts.
- **Flexible**: Allows construction of arbitrary lattice models.
- **Lightweight**: Minimal dependencies (NumPy, SciPy).

## Inputs & Outputs
- **Inputs**:
  - Python scripts defining system parameters (sites, interactions).
  - Explicit Hamiltonian construction calls.
- **Outputs**:
  - NumPy arrays (eigenvalues, eigenvectors).
  - Green's function data.
  - State observables.

## Performance Characteristics
- **Scale**: Limited to small systems due to exponential scaling of ED (typically < 20 sites for full diagonalization).
- **Efficiency**: Uses SciPy sparse matrices for memory efficiency.
- **Parallelization**: Relies on underlying NumPy/SciPy optimization; not explicitly MPI/GPU parallelized.

## Limitations & Known Constraints
- **System Size**: Restricted by the exponential Hilbert space growth of exact diagonalization.
- **Development Status**: Marked as "under development" (alpha/beta stage).
- **Documentation**: Minimal (README and code examples).
- **Stability**: API may change; not a production-grade code for large-scale HPC.

## Comparison with Other Codes
- **vs QuSpin**: QuSpin is more mature, feature-rich, and optimized (C++ backend) for ED. cmpy is simpler and pure Python.
- **vs ALPS**: ALPS provides a broader suite of solvers (DMRG, QMC) and is C++ based.
- **vs PyBinding**: PyBinding is specialized for tight-binding (single particle); cmpy handles many-body interactions.

## Installation
```bash
pip install git+https://github.com/dylanljones/cmpy.git
```
Or verify source matches requirements.

## Verification & Sources
**Primary sources**:
1.  **Repository**: [dylanljones/cmpy on GitHub](https://github.com/dylanljones/cmpy) (Verified ownership and content).
2.  **Topics**: Labeled with `exact-diagonalization`, `hubbard-model`, `greens-functions`.

**Verification status**: ✅ VERIFIED
- **Authenticity**: Confirmed existing repository with relevant physics code.
- **Active**: Last commits within recent history (checked via GitHub interface).
- **Scope**: Small research/utility library, not a major community package.
