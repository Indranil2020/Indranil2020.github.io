# Many-Body-Localization-Exact-Diagonalization (MBL_ED)

## Official Resources
- Homepage: https://github.com/Tcm0/Many-Body-Localization-Exact-Diagonalization
- Source Repository: https://github.com/Tcm0/Many-Body-Localization-Exact-Diagonalization
- License: MIT License

## Overview
**Many-Body-Localization-Exact-Diagonalization** (referred to here as **MBL_ED**) is a specialized Python codebase for studying the **Many-Body Localization (MBL)** transition in 1D spin chains, specifically the random-field Heisenberg model. It utilizes Exact Diagonalization (ED) to compute diagnostics of the MBL transition, such as level statistics (gap ratios) and entanglement entropy. It is designed to efficiently calculate highly excited eigenstates, which are crucial for probing the infinite-temperature physics relevant to MBL.

**Scientific domain**: Many-Body Localization, Quantum Chaos, Disordered Systems
**Target user community**: Researchers investigating ergodicity breaking and thermalization

## Theoretical Methods
- **Shift-Invert Exact Diagonalization**: Uses the shift-invert spectral transformation $(H - \sigma I)^{-1}$ to converge to eigenvalues near a target energy $\sigma$ (often center of spectrum).
- **Level Statistics**: Computation of the adjacent gap ratio $r$ to distinguish between Wigner-Dyson (ergodic) and Poissonian (localized) statistics.
- **Entanglement Entropy**: Calculation of bipartite von Neumann entropy to detect Area-law vs Volume-law scaling.

## Capabilities (CRITICAL)
- **Targeted Solvers**: Optimized for finding interior eigenvalues of sparse matrices, essential for MBL studies.
- **Disorder Averaging**: Scripts set up to run over many disorder realizations.
- **Diagnostics**: Built-in functions for:
  - Level spacing ratios ($r$-parameter).
  - Half-chain Entanglement Entropy.
  - Optical Conductivity (in some versions).

## Key Features

### Specialized MBL Tools:
- **E_inf_temp**: Routines to target energy densities corresponding to infinite temperature states.
- **Disorder Management**: Efficient handling of random field configurations (Box distribution, etc.).

### Python/SciPy Stack:
- **Implementation**: Uses `scipy.sparse.linalg.eigsh` with specific parameters for shift-invert mode.
- **Simplicity**: Concise, readable code ideal for reproducing standard MBL results.

## Inputs & Outputs
- **Input formats**:
  - Python scripts defining system size $L$, disorder strength $W$, and number of realizations.
- **Output data types**:
  - Text/NumPy files containing averaged gap ratios and entropies vs disorder $W$.

## Interfaces & Ecosystem
- **Dependencies**: NumPy, SciPy.
- **Relation**: Similar to functionalities found in `QuSpin`, but as a lightweight, focused repository.

## Workflow and Usage
User modifies parameters in the main script (e.g., `L=12`, `W_list = [...]`) and runs the script. The code loops over disorder realizations, constructs the sparse Hamiltonian, diagonalizes, and accumulates statistics.

## Performance Characteristics
- **System Size**: Typically limited to $L \approx 16-18$ spins for shift-invert ED due to matrix inversion cost.
- **Throughput**: Good for sweeping phase diagrams on cluster nodes.

## Comparison with Other Codes
| Feature | MBL_ED | QuSpin | xdiag |
| :--- | :--- | :--- | :--- |
| **Focus** | MBL / Disorder Averaging | General Many-Body (Dynamics) | General Many-Body (Ground State) |
| **Solver** | Shift-Invert (Interior Eigs) | General (Full/Sparse) | General (Lanczos/Full) |
| **Efficiency** | Optimized for Interior States | General Purpose | Ground State Optimized |
| **Complexity** | Minimal (Script-like) | Moderate (Library) | High (Library) |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/Tcm0/Many-Body-Localization-Exact-Diagonalization

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Niche: focused research code for MBL diagnostics.
