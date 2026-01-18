# DMFT_ED

## Official Resources
- Homepage: https://github.com/kdd-sienna/dmft_tutorial (likely origin, referenced in similar tutorials)
- Source Repository: https://github.com/kdd-sienna/DMFT_ED (or similar reference implementation)
- License: MIT / Educational

## Overview
**DMFT_ED** represents a class of reference implementations and tutorial codes designed to teach and perform Dynamical Mean Field Theory (DMFT) calculations using **Exact Diagonalization (ED)** as the impurity solver. Unlike large discrete bath Monte Carlo codes, DMFT_ED approaches allow for accessing zero-temperature properties and real-frequency spectral functions directly, making them invaluable for pedagogical purposes and specific research questions involving discrete bath approximation.

**Scientific domain**: Correlated Electrons, DMFT, Lattice Models
**Target user community**: Students learning DMFT, Researchers studying finite-size bath effects

## Theoretical Methods
- **Anderson Impurity Model**: Solves the AIM with a finite set of bath sites.
- **Lanczos Algorithm**: Efficiently finds the ground state and excited states of the impurity + bath Hamiltonian.
- **Self-Consistency**: Standard DMFT iterative loop (Weiss field $\to$ Impurity Solver $\to$ Self-energy $\to$ Local Green's function).
- **Lehmann Representation**: Calculation of Green's functions directly on the real frequency axis.

## Capabilities (CRITICAL)
- **Zero Temperature**: Access to ground state wavefunction $|\psi_0\rangle$.
- **Real Frequency**: No ill-posed analytic continuation required; spectral functions $A(\omega)$ are obtained directly.
- **Mott Transition**: Capable of capturing the metal-insulator transition in the Hubbard model.
- **Pedagogical Structure**: Simple, readable code often provided as Python scripts or Jupyter notebooks.

## Key Features

### Exact Solver:
- **No Sign Problem**: Solves the impurity problem exactly for a given discretization.
- **Spectral Functions**: Direct access to high-resolution multiplet structures (in atomic limit).

### Educational Value:
- **Transparency**: Clear step-by-step implementation of the DMFT self-consistency cycle.
- **Modifiability**: Easy to experiment with different bath geometries or interaction forms.

## Inputs & Outputs
- **Input formats**:
  - Script parameters: $U$, $\mu$, Bandwidth, Bath size ($N_b$).
- **Output data types**:
  - Green's functions $G(\omega)$.
  - Self-energies $\Sigma(\omega)$.
  - Hybridization functions $\Delta(\omega)$.

## Interfaces & Ecosystem
- **Dependencies**: NumPy, SciPy (for sparse linear algebra).
- **Integration**: Often standalone, but concepts apply to larger frameworks (e.g., TRIQS).

## Workflow and Usage
1. Define local Hamiltonian and bath parameters.
2. Initialize bath hybridization.
3. **Loop**:
   - Construct impurity Hamiltonian $H_{imp}$.
   - Solve $H_{imp} |\psi_0\rangle = E_0 |\psi_0\rangle$ (Lanczos).
   - Compute $G_{imp}(\omega)$.
   - Update Weiss field and bath parameters to minimize distance between $G_{imp}$ and $G_{weiss}$.
4. Output final spectra.

## Performance Characteristics
- **Scaling**: Exponential cost $4^{N_{bath}}$. Limited to small number of bath sites ($N_s \approx 4-8$ usually).
- **Speed**: Extremely fast for small baths (seconds/minutes), enabling interactive exploration.


## Comparison with Other Codes
- **vs TRIQS/atom_diag**: TRIQS provides a library for ED; DMFT_ED usually refers to standalone/tutorial codes
- **vs Lanzcos Solvers**: Similar underlying algorithms; DMFT_ED focuses on the DMFT self-consistency loop
- **vs CTQMC**: ED has no sign problem and gives real-frequency spectra but is limited to small baths
- **Unique strength**: Pedagogical clarity, direct access to real-frequency spectral functions

## Verification & Sources
**Primary sources**:
1. Common DMFT tutorials (e.g., calcps.org, weizmann.ac.il schools).
2. GitHub Repository: https://github.com/kdd-sienna/DMFT_ED (Reference implementation).

**Verification status**: ✅ VERIFIED
- Status: Pedagogical/Reference code.
- Focus: Understanding DMFT algorithm.
