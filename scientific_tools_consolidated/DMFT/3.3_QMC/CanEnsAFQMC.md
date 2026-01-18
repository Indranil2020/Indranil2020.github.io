# CanEnsAFQMC

## Official Resources
- Homepage: https://github.com/TongSericus/CanEnsAFQMC
- Source Repository: https://github.com/TongSericus/CanEnsAFQMC
- License: MIT License

## Overview
**CanEnsAFQMC** is a Julia-based implementation of the Auxiliary-Field Quantum Monte Carlo (AFQMC) method in the **Canonical Ensemble**. Unlike standard Grand Canonical AFQMC which fixes the chemical potential (and thus only average particle number), this code performs simulations at a strictly fixed particle number $N$. This is critical for small systems, cold atoms in traps, or any scenario where particle number fluctuations are unphysical or undesirable.

**Scientific domain**: Cold Atoms, Quantum Dots, Finite Lattice Systems
**Target user community**: Researchers studying finite quantum systems and Hubbard models

## Theoretical Methods
- **Canonical Ensemble AFQMC**: Uses a recursive projection method (or Fourier projection) to fix total particle number $N$ exactly.
- **Hubbard Model**: Specialized for 2D square lattice Hubbard models.
- **Finite Temperature**: Performs finite-T simulations without chemical potential tuning.

## Capabilities (CRITICAL)
- **Fixed Particle Number**: Simulations are strictly at fixed $N$, eliminating the need to search for chemical potential $\mu$.
- **Convergence**: Often converges to ground state properties faster than Grand Canonical methods for finite systems.
- **Observables**: Measures Energy, Momentum Distribution $n(k)$, Spin Correlations, and Entanglement Entropy (Renyi/Von Neumann).
- **Lattices**: Square lattices (1D chains, 2D planes).

## Key Features

### Julia Implementation:
- **Modern Codebase**: Written in pure Julia for readability and performance.
- **Parallelization**: Supports multi-threading/processing via Julia's native parallelism.

### Unique Algorithms:
- **Recursive Projection**: Efficient algorithm to project onto fixed $N$ subspace during propagation.

## Inputs & Outputs
- **Input formats**:
  - Julia parameter scripts (Lx, Ly, U, Beta, n_particles).
- **Output data types**:
  - Text/HDF5 files with measured observables.

## Interfaces & Ecosystem
- **Dependencies**: Julia `LinearAlgebra`, `Statistics`.
- **Standalone**: operates independently but results can be compared with `MonteCarlo.jl` (Grand Canonical).

## Workflow and Usage
```julia
using CanEnsAFQMC
# Define parameters
p = Parameters(Lx=4, Ly=4, U=4.0, beta=10.0, N=14)
# Run simulation
run_afqmc(p)
```

## Performance Characteristics
- **Speed**: $O(N^3)$ or $O(N^4)$ depending on projection method scaling vs Grand Canonical.
- **Accuracy**: Removes systematic errors associated with $\mu$ tuning and particle fluctuation in finite systems.

## Comparison with Other Codes
| Feature | CanEnsAFQMC | MonteCarlo.jl | ALF |
| :--- | :--- | :--- | :--- |
| **Ensemble** | Canonical (Fixed N) | Grand Canonical (Fixed $\mu$) | Grand Canonical (Fixed $\mu$) |
| **Language** | Julia | Julia | Fortran 2003 |
| **Primary Model** | Hubbard (Finite) | Classic/Quantum Spin, Hubbard | Lattice Fermions (General) |
| **Key Advantage** | Exact fixed particle number, no $\mu$-tuning | Versatility, multi-model | Production stability, finite T |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/TongSericus/CanEnsAFQMC
2. Publications by the authors on Canonical Ensemble AFQMC methods (e.g., in PRB or PRE).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Methodology: Canonical projection is a well-established but specialized technique.
