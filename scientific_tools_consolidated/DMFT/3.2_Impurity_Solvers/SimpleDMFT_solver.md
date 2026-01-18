# SimpleDMFT_solver (romainfd)

## Official Resources
- Source Repository: https://github.com/romainfd/DMFT_solver
- License: No explicit license found (Copyright retained by author)
- Language: Python, Jupyter Notebook

## Overview
This repository contains a simple, educational implementation of a DMFT solver. It focuses on the Iterated Perturbation Theory (IPT) method for the single-orbital Hubbard model on a Bethe lattice. It serves as an accessible entry point for understanding the structure and implementation of a DMFT self-consistency loop and perturbative solvers, utilizing a clear Jupyter Notebook format.

**Scientific domain**: Educational Physics, Many-Body Theory
**Target user community**: Students, Beginners in DMFT

## Theoretical Methods
- Dynamical Mean-Field Theory (DMFT)
- Iterated Perturbation Theory (IPT)
- Hubbard Model (Half-filling focus)
- Bethe Lattice Density of States

## Capabilities
- **Exact Solution (Limit)**: Solves the single-band Bethe lattice Hubbard model.
- **Self-Consistency**: Demonstrates the full DMFT loop implementation in Python.
- **IPT Solver**: Implements the efficient IPT approximation for the impurity problem.
- **Spectral Functions**: Calculates Green's functions and self-energies on real frequencies.

## Key Strengths
### Simplicity:
- Minimal codebase, easy to read and modify.
### Educational Value:
- Demonstrates the core logic of DMFT without the complexity of optimized production codes.
- Jupyter Notebook format allows for interactive learning and plotting.

## Inputs & Outputs
- **Input parameters**:
  - Model parameters defined in notebook: $U$ (interaction), $\beta$ (inverse temperature), $D$ (bandwidth).
  - Code parameters: `n_loops`, `mixing`.
- **Outputs**:
  - Plots of Green's functions and Self-energies vs frequency.
  - Quasiparticle weight $Z$.

## Interfaces & Ecosystem
- **Language**: Python (Primary), Jupyter Notebook.
- **Dependencies**: NumPy, Matplotlib, SciPy.

## Performance Characteristics
- **Speed**: Very fast (IPT is analytical/algebraic).
- **Compute**: Runs instantly on standard laptops.

## Limitations & Known Constraints
- **Scope**: Limited to single orbital, Bethe lattice (infinite coordination number).
- **Solver**: IPT is approximate (strictly valid at half-filling for this implementation).
- **Production Use**: Not intended for material science production runs.

## Comparison with Other Codes
- **vs TRIQS/w2dynamics**: This is a toy code for learning, not a library.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/romainfd/DMFT_solver
2. Description: "LISA DMFT solver"

**Verification status**: ✅ VERIFIED
- Source code: OPEN
- Purpose: Educational / "Toy Code"
