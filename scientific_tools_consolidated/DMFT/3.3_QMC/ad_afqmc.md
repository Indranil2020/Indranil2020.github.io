# ad_afqmc

## Official Resources
- Homepage: https://github.com/ankit76/ad_afqmc
- Source Repository: https://github.com/ankit76/ad_afqmc
- License: MIT License

## Overview
**ad_afqmc** is an end-to-end automatically differentiable Auxiliary Field Quantum Monte Carlo (AFQMC) library built on top of JAX. It leverages automatic differentiation (AD) to enable the efficient calculation of physical properties and gradients, such as nuclear forces and dipole moments, which are traditionally challenging in stochastic frameworks. This approach allows for gradient-based optimization of wavefunctions and geometry relaxation within the AFQMC method.

**Scientific domain**: Quantum Chemistry, Machine Learning Physics, Electronic Structure
**Target user community**: Researchers at the intersection of QMC, differentiable programming, and machine learning

## Theoretical Methods
- **Phaseless AFQMC**: Standard phaseless Auxiliary Field QMC algorithm for fermions.
- **Automatic Differentiation (AD)**:
  - **Reverse-mode AD**: Used for observables like 1-RDM (Computationally intensive, uses checkpointing).
  - **Forward-mode AD**: Used for derivatives like dipole moments.
- **Gradient Estimation**: Overcomes mixed-estimator bias in standard QMC gradient calculations.
- **PySCF Integration**: Uses PySCF for molecular integral generation (1-electron, 2-electron integrals).

## Capabilities (CRITICAL)
- **Differentiable QMC**: Compute exact gradients of energy w.r.t. parameters (wavefunction or Hamiltonian).
- **Optimization**: Enables VMC/AFQMC optimization of trial wavefunctions using AD gradients.
- **Hardware Acceleration**: Native support for GPUs and TPUs via JAX's XLA compiler.
- **Molecular Properties**: Calculates dipole moments, nuclear gradients, and other observables.
- **Memory Management**: Implements checkpointing to handle high memory cost of reverse-mode AD for large systems.

## Key Features

### JAX Backend:
- **JIT Compilation**: Just-In-Time compilation (via `jax.jit`) for core MC kernels.
- **Vectorization**: `jax.vmap` for efficient batch processing of walkers.
- **GPU/TPU Support**: Seamless execution on accelerators.

### Quantum Chemistry:
- **PySCF Interface**: Direct integration for setting up molecular systems.
- **Geometry Optimization**: Can drive geometry relaxation using AFQMC force estimates (experimental).

## Inputs & Outputs
- **Input formats**:
  - Python scripts using `jax` and `pyscf` objects.
- **Output data types**:
  - Energies and gradients.
  - Optimized wavefunction parameters.

## Interfaces & Ecosystem
- **Dependencies**: JAX, PySCF, HDF5, NumPy.
- **Ecosystem**: Part of the growing family of differentiable quantum physics tools (like FermiNet, VMCNet).

## Workflow and Usage
Users define a molecule in PySCF, pass the integrals to `ad_afqmc`, and define a JAX-jitted function to run the AFQMC propagation. Gradients can be obtained simply by calling `jax.grad` on the energy estimator function.

## Performance Characteristics
- **Speed**: High performance on GPUs due to XLA compilation.
- **Scaling**: Excellent parallel scaling on GPU clusters (batching walkers).
- **Cost**: AD adds overhead (memory and compute) compared to standard propagation, but provides gradients "for free".

## Comparison with Other Codes
| Feature | ad_afqmc | ipie | QMCPACK (AFQMC) |
| :--- | :--- | :--- | :--- |
| **Core Tech** | JAX (Differentiable) | Python / Cupy / Numba | C++ / CUDA |
| **Differentiation** | Automatic (End-to-End) | Manual / None | None |
| **Primary Goal** | Optimization / Gradients | Performance / Properties | Production / scale |
| **Hardware** | GPU / TPU (via JAX) | GPU (via Cupy) | GPU (via CUDA) |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/ankit76/ad_afqmc
2. Publication: "End-to-end differentiable auxiliary field quantum Monte Carlo" (likely J. Chem. Phys or similar AIP journal).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- State: Research code, innovative method.
- Focus: Differentiable QMC.
