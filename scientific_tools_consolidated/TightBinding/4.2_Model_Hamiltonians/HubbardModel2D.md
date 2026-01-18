# HubbardModel2D

## Official Resources
- **Repository**: https://github.com/ryanlevy/HubbardModel2D
- **License**: MIT License

## Overview
**HubbardModel2D** is a specialized C++ code for performing **Exact Diagonalization (ED)** of the single-orbital Fermi-Hubbard model on 1D chains and 2D square lattices. It uses the **Lanczos algorithm** to compute ground state energies, wavefunctions, and spectral properties. While primarily designed for small clusters due to the exponential scaling of ED, it serves as a transparent and efficient tool for benchmarking and studying strong correlation effects on finite lattices.

**Scientific domain**: Strongly Correlated Electrons, Quantum Magnetism
**Target user community**: Students and researchers characterizing the Hubbard model

## Theoretical Methods
- **Fermi-Hubbard Model**: $H = -t \sum_{\langle ij \rangle \sigma} (c_{i\sigma}^\dagger c_{j\sigma} + h.c.) + U \sum_i n_{i\uparrow} n_{i\downarrow}$.
- **Exact Diagonalization**: Construction of the full sparse Hamiltonian in the occupation number basis for fixed particle number and spin.
- **Lanczos Algorithm**: Iterative method to find the lowest eigenvalues and eigenvectors (Ground State).
- **Conserved Quantities**: Exploits $N_{\uparrow}$ and $N_{\downarrow}$ conservation (but not spatial symmetries).

## Capabilities
- **Simulations**:
  - Ground State Energy $E_0$.
  - Wavefunction amplitues.
  - Spectral Functions $A(\omega)$ (via Spectra library).
- **Geometries**:
  - 1D Chains (PBC/OBC).
  - 2D Square Clusters ($2 \times 2$, $2 \times 4$, $4 \times 4$, etc.).
- **Observables**:
  - Spin-spin correlations $\langle S_i \cdot S_j \rangle$.
  - Charge correlations $\langle n_i n_j \rangle$.

## Key Strengths
- **Simplicity**: Codebase is relatively small and uses the highly readable `Eigen` library for linear algebra, making it easy to modify.
- **Backend Options**: Can switch between `ietl` (classic) and `Spectra` (modern) for the Lanczos implementation.
- **Benchmarking**: Ideal for generating exact reference data for verifying approximate methods like VMC or DMRG on small systems.

## Inputs & Outputs
- **Inputs**: C++ main file modification (recompilation required for changes) to set $t$, $U$, lattice size.
- **Outputs**:
  - Console output of energies.
  - Text files for correlation functions.

## Interfaces & Ecosystem
- **Dependencies**: Eigen3, Spectra (header-only).
- **Ecosystem**: Standalone C++ tool.

## Performance Characteristics
- **Scaling**: Exponential Hilbert space growth. Feasible for up to $\sim 16$ sites (e.g., $4 \times 4$ at half-filling is pushing memory limits).
- **Parallelism**: OpenMP shared-memory parallelism for matrix-vector multiplication.

## Comparison with Other Codes
- **vs. ALPS**: ALPS handles full spatial symmetries, extending the reach to slightly larger systems ($\sim 20$ sites). HubbardModel2D is simpler but restricted to smaller sizes.
- **vs. EDLib**: Similar scope; HubbardModel2D is an application for the Hubbard model, while EDLib is a library.

## Application Areas
- **Mott Physics**: Studying the metal-insulator transition in small clusters.
- **Pairing Correlations**: Checking for d-wave pairing signatures in the 2D Hubbard ground state.

## Community and Support
- **Development**: Ryan Levy.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/ryanlevy/HubbardModel2D](https://github.com/ryanlevy/HubbardModel2D)
- **Verification status**: ✅ VERIFIED
  - Functional educational/research code.
