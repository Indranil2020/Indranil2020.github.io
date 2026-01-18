# Heisenberg

## Official Resources
- **Repository**: https://github.com/muammar/heisenberg
- **License**: MIT License

## Overview
**Heisenberg** is a lightweight Python program designed for the exact analysis of quantum spin chains. It constructs the Hamiltonian matrix for the **Spin-1/2 Heisenberg model** in the $S^z$ basis, diagonalizes it to extract the full spectrum (eigenvalues and eigenvectors), and computes advanced quantities like the **Total Position Spread (TPS)** tensor, which is relevant for studying localization and insulating behavior.

**Scientific domain**: Quantum Magnetism, Educational Physics
**Target user community**: Students and researchers studying small 1D spin systems

## Theoretical Methods
- **Exact Diagonalization (ED)**: Full diagonalization of the Hamiltonian matrix $H$.
- **Heisenberg Model**: $H = \sum_{i} J \mathbf{S}_i \cdot \mathbf{S}_{i+1}$ (supports Antiferromagnetic/Ferromagnetic couplings).
- **Boundary Conditions**:
  - Open Boundary Conditions (OBC).
  - Periodic Boundary Conditions (PBC).
- **Position Spread**: Calculation of the second moment of the position operator (TPS tensor) to characterize charge/spin fluctuation.

## Capabilities
- **Hamiltonian**: Generation of sparse matrices for $N$ spins.
- **Spectrum**: Complete set of energy levels and eigenstates.
- **Observables**:
  - Ground state energy.
  - TPS tensor (for distinguishing metals/insulators in electronic mappings).
- **Simplicity**: Pure Python implementation using SciPy sparse matrices.

## Key Strengths
- **Accessibility**: Less than 500 lines of clear Python code, making it an excellent learning resource for understanding how ED works "under the hood."
- **TPS Tensor**: One of the few simple codes that explicitly implements the Total Position Spread calculation out-of-the-box.

## Inputs & Outputs
- **Inputs**:
  - System size $N$.
  - Coupling constant $J$.
  - Boundary condition flags.
- **Outputs**:
  - Printed eigenvalues.
  - Plots of eigenvalues and TPS scaling.

## Interfaces & Ecosystem
- **Dependencies**: NumPy, SciPy, Matplotlib.
- **Format**: Standalone script.

## Performance Characteristics
- **Scaling**: Exponential $2^N$. Practical limit on a laptop is $N \approx 14-16$ spins.
- **Memory**: Stores the Hamiltonian, limiting size compared to Lanczos-based codes (like EDLib) that don't need the full matrix.

## Comparison with Other Codes
- **vs. QuSpin**: QuSpin is a professional research framework for ED; Heisenberg is a minimal script. Use QuSpin for serious research; use Heisenberg for learning or quick checks.
- **vs. JHeisenbergED**: Similar scope, but written in Python instead of Julia.

## Application Areas
- **Pedagogy**: Teaching quantum mechanics of many-body systems.
- **Localization**: Testing criteria for localization via the TPS tensor.

## Community and Support
- **Development**: Muammar El Khatib.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/muammar/heisenberg](https://github.com/muammar/heisenberg)
- **Verification status**: ✅ VERIFIED
  - Functional Python script for small-scale ED.
