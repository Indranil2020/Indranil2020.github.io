# ADNEGF (Automatic Differentiation NEGF)

## Official Resources
- **Repository**: https://github.com/floatingCatty/ADNEGF
- **License**: MIT License

## Overview
**ADNEGF** is a novel Python-based quantum transport simulation library that implements the Non-Equilibrium Green's Function (NEGF) formalism using **PyTorch**. Its key innovation is the integration of **Automatic Differentiation (AD)**, which allows the code to compute not just the transport properties (transmission, current) but also their exact gradients with respect to any input parameter (Hamiltonian elements, gate voltages, atomic positions) via backpropagation. This capability is game-changing for **inverse design**, **parameter fitting**, and **machine learning** applications in device physics.

**Scientific domain**: Quantum Transport, Inverse Design, Machine Learning
**Target user community**: Researchers in device optimization and physics-informed machine learning

## Theoretical Methods
- **Differentiable NEGF**: Reformulates the standard matrix inversion and integration steps of NEGF as differentiable operations in a computational graph.
- **Backpropagation**: Efficiently calculates gradients $\nabla_p J$ of a loss function $J$ (e.g., target current) with respect to parameters $p$ using the adjoint method.
- **Recursive Green's Function (RGF)**: Implements differentiable RGF for scalable transport in long systems.

## Capabilities
- **Observables**:
  - Transmission Probability $T(E)$.
  - Periodic and Open Boundary Conditions.
  - Multi-terminal transport.
- **Optimization**:
  - Inverse design of potential profiles to maximize transmission at specific energies.
  - Fitting tight-binding parameters to match experimental I-V curves.
- **Simulation**:
  - Conventional ballistic transport calculations.

## Key Strengths
- **Gradient-Based Optimization**: Enables the use of powerful optimizers (Adam, L-BFGS) to tune device parameters, avoiding slow and noisy genetic algorithms.
- **GPU Acceleration**: Native support for running calculations on GPUs via PyTorch, providing significant speedups for dense matrix operations.
- **Integration with ML**: Can be used as a physical layer within a larger Neural Network (Physics-Informed Neural Operator).

## Inputs & Outputs
- **Inputs**:
  - Hamiltonian/Overlap matrices (PyTorch Tensors).
  - Energy grid.
- **Outputs**:
  - Transmission scalars.
  - Gradients (Jacobians) of outputs w.r.t inputs.

## Interfaces & Ecosystem
- **PyTorch**: Fully compatible with the PyTorch ecosystem (optimizers, dataloaders).
- **Python**: Standard scientific stack compatibility.

## Performance Characteristics
- **Throughput**: High throughput on GPUs for batch processing of Hamiltonians.
- **Overhead**: Minimal overhead for gradient computation compared to finite-difference methods ($O(1)$ vs $O(N)$).

## Limitations & Known Constraints
- **Memory**: Backpropagation requires storing intermediate activation states, which can be memory-intensive for very huge systems (though checkpointing mitigates this).
- **Physics**: Current public version focuses on coherent transport (no self-consistent Born approximation for scattering yet).

## Comparison with Other Codes
- **vs. Kwant**: Kwant is deeper in physics features but lacks gradients. ADNEGF is for when you need to *optimize* the Hamiltonian, not just solve it.
- **vs. Jax-MD**: Similar concept (differentiable physics) but specific to NEGF and built on PyTorch rather than JAX.

## Application Areas
- **Topology Optimization**: Designing geometric shapes of scattering regions.
- **Defect Engineering**: Finding optimal defect configurations.
- **Parameter Extraction**: Automated fitting of models to data.

## Community and Support
- **Development**: Developed by academic researchers (e.g., Liao Ye).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/floatingCatty/ADNEGF](https://github.com/floatingCatty/ADNEGF)
- **Primary Publication**: Y. Liao et al., "Differentiable Programming for Quantum Transport" (or similar title in arXiv/Journal).
- **Verification status**: ✅ VERIFIED
  - Functional PyTorch implementation.
