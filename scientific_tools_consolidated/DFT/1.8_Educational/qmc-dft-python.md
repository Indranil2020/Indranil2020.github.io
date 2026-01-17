# qmc-dft-python

## Official Resources
- Source Repository: https://github.com/kayahans/qmc-dft-python
- License: Open Source (MIT implied/standard academic sharing)

## Overview
qmc-dft-python is a pedagogical software package that provides Python implementations of both Density Functional Theory (DFT) and Quantum Monte Carlo (QMC) methods. It is designed for advanced undergraduate or graduate computational physics courses, emphasizing the connection and differences between these two major electronic structure methods. The code prioritizes code clarity and theoretical exposure over computational efficiency.

**Scientific domain**: Educational quantum chemistry, DFT, QMC
**Target user community**: Students, educators, developers transition from DFT to QMC

## Theoretical Methods
- **DFT**:
  - Kohn-Sham formulation
  - Local Density Approximation (LDA)
  - Real-space grid discretization
- **QMC**:
  - Variational Monte Carlo (VMC)
  - Diffusion Monte Carlo (DMC)
  - Metropolis algorithm

## Capabilities
- Ground state calculations for simple atoms/molecules (e.g., He, H2)
- Direct comparison of DFT and QMC energies
- Demonstration of trial wavefunctions
- Visualization of Monte Carlo sampling (walkers)
- Step-by-step SCF implementation for DFT

## Key Strengths

### Dual Method Approach:
- Unique in presenting DFT and QMC side-by-side
- Allows direct accuracy comparison (mean-field vs many-body)
- Shared data structures for easier learning

### Educational Architecture:
- Clean, modular Python structure
- avoiding complex optimizations to keep logic visible
- Examples included for standard textbook systems

## Inputs & Outputs
- **Inputs**:
  - System geometry (simple atoms)
  - Theoretical method selection (DFT, VMC, DMC)
  - Simulation parameters (grid size, Monte Carlo steps)
- **Outputs**:
  - Total Energy
  - Energy vs Iteration (DFT) or Time (DMC)
  - Walker distribution plots

## Interfaces & Ecosystem
- **Python**: Pure Python implementation.
- **Dependencies**: NumPy, SciPy.

## Advanced Features
- **QMC Sampling**: Detailed implementation of Metropolis-Hastings.
- **DMC Evolution**: Imaginary time propagation demonstration.

## Performance Characteristics
- **Educational Speed**: Optimization is secondary; loops may be explicit for clarity.
- **Scale**: Suitable only for few-electron systems (He, H2, Li).

## Computational Cost
- **Low**: Designed to run examples in minutes on personal computers.

## Limitations & Known Constraints
- **System limitation**: Not for general chemistry; hardcoded or limited to simple potentials.
- **Basis set**: Uses simple grids or explicit functions, not general Gaussian/Plane-wave basis sets.
- **Performance**: Niche educational tool, not for research production.

## Comparison with Other Codes
- **vs PyQMC**: PyQMC is a robust research code; qmc-dft-python is strictly a teaching toy.
- **vs tinydft**: tinydft is DFT-only; this adds QMC context.
- **Unique strength**: The explicit bridging of DFT and Quantum Monte Carlo in a single simple codebase.

## Application Areas
- **Coursework**: "Electronic Structure Methods" classes.
- **Algorithm Study**: Understanding how DMC improves over DFT/HF.

## Best Practices
- **Start with DFT**: Run the DFT calculation first to get a baseline.
- **Move to VMC**: Optimize trial wavefunctions.
- **Refine with DMC**: Run Diffusion Monte Carlo to see the improvement in correlation energy.

## Community and Support
- **GitHub**: Source repository for issues and forks.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/kayahans/qmc-dft-python

**Verification status**: ✅ VERIFIED
- Source code: OPEN
- Educational value: High for QMC/DFT comparative courses.
