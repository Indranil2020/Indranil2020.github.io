# exactdiag (mikeschmitt)

## Official Resources
- Homepage: https://github.com/mikeschmitt/exactdiag
- Source Repository: https://github.com/mikeschmitt/exactdiag
- License: MIT License

## Overview
**exactdiag** is a Python package for performing exact diagonalization of fermionic many-body systems. Uniquely, it leverages **Numba** to Just-In-Time (JIT) compile critical inner loops, allowing it to achieve performance close to compiled languages like C++ or Fortran while efficient Python scripting. It focuses on the Anderson Impurity Model and Hubbard models, providing tools for Green's functions and spectral densities.

**Scientific domain**: Correlated Electrons, DMFT Impurity Solvers
**Target user community**: Python users requiring fast ED for small fermionic clusters

## Theoretical Methods
- **Fermionic ED**: Creation/Annihilation operator algebra treating fermion signs correctly.
- **Lehmann Representation**: Calculation of dynamical quantities ($G(\omega)$) via spectral expansion.
- **Linear Operators**: Uses `scipy.sparse.linalg.LinearOperator` to avoid storing full matrices when possible.

## Capabilities (CRITICAL)
- **Numba Acceleration**: JIT compilation of Hamiltonian action on state vectors.
- **Model Support**: Specialized for Single-Impurity Anderson Models (SIAM) and Hubbard chains/clusters.
- **Spectral Functions**: Efficient computation of Zero-temperature Green's functions.
- **Symmetries**: Can exploit particle number and spin conservation.

## Key Features

### Performance:
- **Python + Numba**: The ease of Python with the speed of machine code.
- **Sparse Algebra**: Interface with generic SciPy eigensolvers (ARPACK).

### Usability:
- **Object-Oriented**: Classes for `Basis`, `Hamiltonian`, and `GreenFunction`.

## Inputs & Outputs
- **Input formats**:
  - Python scripts defining hopping $t_{ij}$ and interaction $U_{ijkl}$.
- **Output data types**:
  - Arrays of frequencies and spectral weights.

## Interfaces & Ecosystem
- **Dependencies**: Python, NumPy, SciPy, Numba.
- **Integration**: Can be used as a lightweight solver in a Python-based DMFT loop.

## Workflow and Usage
```python
import exactdiag as ed
# Define basis
basis = ed.Basis(nsites=4, nup=2, ndn=2)
# Create Hamiltonian
H = ed.Hamiltonian(basis, ...parameters...)
# Diagnose
evals, evecs = H.diagonalize()
```

## Performance Characteristics
- **Speed**: Significantly faster than pure Python/NumPy implementations due to Numba.
- **Memory**: Standard ED limits apply, but sparse formulation helps.

## Comparison with Other Codes
| Feature | exactdiag (mikeschmitt) | xdiag (awietek) | TRIQS |
| :--- | :--- | :--- | :--- |
| **Language** | Python + Numba | C++ / Julia | C++ / Python |
| **Focus** | Impurity Solvers (SIAM) | General Lattice ED | DMFT / Green's Functions |
| **Speed** | Numba JIT (Fast) | C++ Optimized (Very Fast) | C++ Optimized (Very Fast) |
| **Usage** | Lightweight Python | HPC / Production | Heavy Framework |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/mikeschmitt/exactdiag

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Focus: Numba-accelerated ED.
