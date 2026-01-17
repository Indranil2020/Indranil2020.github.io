# GradDFT

## Official Resources
- Homepage: https://github.com/XanaduAI/GradDFT
- Source Repository: https://github.com/XanaduAI/GradDFT
- License: Apache License 2.0

## Overview
GradDFT (Gradient-based Density Functional Theory) is a pioneering JAX-based library for Machine Learning enhanced DFT (ML-DFT). It leverages the power of Differentiable Programming (Diff-Prog) to enable the training of Exchange-Correlation (XC) functionals and basis sets using gradient descent. By backpropagating through the entire Self-Consistent Field (SCF) cycle, GradDFT allows researchers to optimize neural network functionals directly against high-accuracy data (like CCSD(T) densities or energies), opening new frontiers in functional design.

**Scientific domain**: ML-DFT, Differentiable Programming, Quantum Chemistry
**Target user community**: Functional developers, AI/ML researchers in chemistry

## Theoretical Methods
- **Kohn-Sham DFT**: Standard formulation.
- **Differentiable Programming**: AD (Automatic Differentiation) through the SCF loop.
- **Neural XC Functionals**: Representing $E_{xc}[\rho]$ with neural networks.
- **JAX**: High-performance numerical computing with composable transformations.
- **Basis Set Optimization**: Differentiating w.r.t Gaussian exponents/coefficients.

## Capabilities (CRITICAL)
- **Trainable Functionals**: Define $E_{xc}$ as a Neural Network and train it.
- **Differentiable SCF**: Obtain derivatives of Energy/Density w.r.t any parameter.
- **Properties**: Energy, Forces, Dipoles, Densities.
- **Hardware Acceleration**: Native support for GPUs and TPUs via JAX.
- **Integration**: Uses PySCF for integrals and initial guesses (in some workflows) or native JAX implementations.

## Key Strengths

### End-to-End Differentiability:
- Unlike "train-then-freeze" approaches, GradDFT allows optimization of parameters by determining how they affect the *final* converged density.
- Enables "Physics-Informed" constraints (e.g., exact limits).

### JAX Ecosystem:
- Just-In-Time (JIT) compilation for extreme speed.
- Vectorization (vmap) for batched calculations of molecules.

### Functional Freedom:
- Users can write arbitrary Python code for the functional, and JAX handles the gradients. No restricted forms.

## Inputs & Outputs
- **Inputs**:
  - Molecular geometry (Atom positions, Types).
  - Neural Network architecture (Flax/Haiku).
  - Training targets (Energies, Densities).
- **Outputs**:
  - Optimized NN weights.
  - Loss histories.
  - Final Functional object.

## Interfaces & Ecosystem
- **JAX**: Core framework.
- **PySCF**: Often used as the integral engine or benchmark comparison.
- **Optax**: For optimization algorithms (Adam, SGD).

## Advanced Features
- **Solid State**: Experimental support for periodic boundary conditions (depending on version).
- **Automatic Basis Optimization**: Refining basis sets for specific chemical environments.

## Performance Characteristics
- **Training Speed**: Highly efficient on GPUs due to JAX XLA compilation.
- **Memory**: Backpropagation through SCF (if unrolled) can be memory intensive; uses implicit differentiation (Deep Equilibrium Models) techniques to mitigate this.

## Computational Cost
- **High (Training)**: Training a functional is costly.
- **Varies (Inference)**: Depends on the NN size; typically comparable to hybrid DFT.

## Limitations & Known Constraints
- **Maturity**: Research code; API may change.
- **Features**: Lacks full feature set of VASP/Gaussian (e.g., extensive post-processing, solvation). It is a *laboratory* for functionals.
- **Complexity**: Requires understanding of JAX and AD concepts.

## Comparison with Other Codes
- **vs DQC**: DQC (Differentiable Quantum Chemistry) is PyTorch-based; GradDFT is JAX-based. JAX often offers superior performance for scientific computing primitives.
- **vs PySCF**: PySCF is the workhorse for standard DFT; GradDFT uses PySCF principles but adds the AD layer for *learning* DFT.
- **vs DeepH**: DeepH learns the Hamiltonian; GradDFT learns the Functional.
- **Unique strength**: The most accessible framework for Differentiable DFT research on the JAX stack.

## Application Areas
- **Functional Design**: Creating "universal" NN functionals.
- **System-Specific DFT**: Training a lightweight functional for a specific protein or polymer.
- **Density Inversion**: Finding the exact KS potential corresponding to an accurate density.

## Best Practices
- **Use GPUs**: JAX on CPU is slow; JAX on GPU is flying.
- **Implicit Differentiation**: Use implicit gradient modes to avoid storing the entire SCF trajectory graph.
- **Start Simple**: Train an LDA-like NN before moving to non-local or orbital-dependent forms.

## Community and Support
- **GitHub**: Managed by Xanadu AI (Quantum Computing company).
- **Documentation**: Tutorials and API docs available.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/XanaduAI/GradDFT
2. Documentation: Xanadu website/docs.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (Apache 2.0)
- Innovation: Leading tool in "Differentiable DFT".
