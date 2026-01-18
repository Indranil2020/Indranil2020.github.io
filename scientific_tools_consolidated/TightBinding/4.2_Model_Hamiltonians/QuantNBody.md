# QuantNBody

## Official Resources
- **Homepage**: https://quantnbody.github.io/QuantNBody/
- **Repository**: https://github.com/SYalouz/QuantNBody
- **License**: GPL-3.0

## Overview
**QuantNBody** is a Python package designed to facilitate the manipulation of many-body operators and wave functions in the second quantization formalism. Developed effectively as a "quantum chemistry playground," it allows users to implement and test various electronic structure methods—such as Configuration Interaction (CI) and Coupled Cluster (CC)—from scratch. It is particularly valuable for educational purposes and for method developers prototyping new many-body algorithms.

**Scientific domain**: Quantum Chemistry, Many-Body Physics
**Target user community**: Educators, students, and developers of new quantum algorithms

## Theoretical Methods
- **Second Quantization**: Representation of operators using creation/annihilation strings.
- **Basis States**: Management of fermionic Slater determinants.
- **Hamiltonian**: 1-body and 2-body integral representation ($h_{pq}$, $g_{pqrs}$).
- **Solvers**:
  - Full Configuration Interaction (FCI).
  - Complete Active Space (CAS-CI).
  - Restricted Active Space (RAS-CI).

## Capabilities
- **Operator Algebra**: Tools to build arbitrary spin-conserving or spin-flipping operators.
- **Matrix Generation**: Automatic construction of the Hamiltonian matrix in a specific many-body basis.
- **Interfaces**: Can read integrals from **Psi4** to perform calculations on real molecules.
- **Analysis**: Calculation of reduced density matrices (1-RDM, 2-RDM).

## Key Strengths
- **Transparency**: Unlike black-box chemistry codes, QuantNBody exposes the raw matrix vector operations, making it excellent for understanding *how* FCI or CASSCF works.
- **Prototyping**: Ideal for testing hybrid quantum-classical algorithms (e.g., VQE ansätze) where access to the explicit Hamiltonian matrix is needed.
- **Python-Native**: Fully integrated with the SciPy stack.

## Inputs & Outputs
- **Inputs**:
  - 1- and 2-electron integrals (NumPy arrays).
  - Active space definitions.
- **Outputs**:
  - Ground and excited state energies.
  - Wavefunction coefficients ($\mathbf{C}$ vector).

## Interfaces & Ecosystem
- **Psi4**: Used as a driver to generate molecular integrals.
- **PySCF**: Can be used similarly for integrals.

## Performance Characteristics
- **Scaling**: Limited by the combinatorial explosion of the Hilbert space (FCI limit). Useful for small active spaces (e.g., 10-14 orbitals).
- **Efficiency**: Uses efficient bit-manipulation and sparse matrices where possible, but not a replacement for production-level codes like ORCA or Molpro.

## Comparison with Other Codes
- **vs. PySCF**: PySCF is a full production suite. QuantNBody is a lightweight library for manipulating the algebraic structure of the problem.
- **vs. Hande-QMC**: Hande is for QMC; QuantNBody is for determinstic CI methods.

## Application Areas
- **Education**: Teaching many-body theory courses.
- **Algorithm Dev**: Testing novel embedding schemes or quantum computing mappings.

## Community and Support
- **Development**: Saad Yalouz.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/SYalouz/QuantNBody](https://github.com/SYalouz/QuantNBody)
- **Primary Publication**: S. Yalouz et al., J. Open Source Softw. 7(80), 4759 (2022).
- **Verification status**: ✅ VERIFIED
  - Active educational tool.
