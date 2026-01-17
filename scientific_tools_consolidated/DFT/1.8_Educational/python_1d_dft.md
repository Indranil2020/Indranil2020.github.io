# python_1d_dft

## Official Resources
- Source Repository: https://github.com/tamuhey/python_1d_dft
- License: MIT License

## Overview
python_1d_dft is a minimalistic, educational density functional theory code implemented in Python. It simulates a 1D harmonic oscillator system to demonstrate the fundamental concepts of DFT, including the Kohn-Sham equations, local density approximation (LDA), and self-consistent field (SCF) cycles. It is designed specifically for students and beginners to understand the internal mechanics of a DFT calculation without the complexity of a full-scale production code.

**Scientific domain**: Educational theory, 1D model systems
**Target user community**: Students, beginners in computational physics, educators

## Theoretical Methods
- Kohn-Sham Density Functional Theory (KS-DFT)
- Local Density Approximation (LDA)
- 1D Harmonic Oscillator potential
- Real-space finite difference discretization
- Self-consistent field (SCF) iteration
- Eigenvalue conceptual demonstration

## Capabilities
- Solves 1D Schrödinger equation (Kohn-Sham)
- Calculates total ground state energy
- Visualizes electron density
- Demonstrates convergence behavior
- Modular Python implementation for easy reading

## Key Strengths

### Educational Clarity:
- <100 lines of core logic
- Heavily commented code
- Focus on readability over performance
- Isolates key DFT steps (Hamiltonian construction, diagonalization, density update)

### Pure Python:
- No compilation required
- Uses standard libraries (NumPy, SciPy, Matplotlib)
- Easy to modify and experiment with

## Inputs & Outputs
- **Input parameters**:
  - Number of electrons
  - Grid points
  - Mixing parameter
- **Output data**:
  - Total energy
  - Eigenvalues
  - Density plots (matplotlib)
  - Convergence history

## Interfaces & Ecosystem
- **Python**: Native Python code, integrates with NumPy/Matplotlib
- **Jupyter**: Well-suited for interactive notebook tutorials

## Advanced Features
- **Visualization**: Built-in plotting for density and potential
- **Simplicity**: Can be easily extended to other 1D potentials by the user

## Performance Characteristics
- **Speed**: Instantaneous for model systems
- **System size**: Limited to simple 1D models
- **Parallelization**: Serial only (educational)

## Computational Cost
- **Minimal**: Runs on any standard laptop or Google Colab instance in seconds.

## Limitations & Known Constraints
- **1D Only**: Restricted to one-dimensional model systems.
- **Model Potentials**: Not for real materials or molecules.
- **Educational**: Not performance, highly unoptimized for large grids.

## Comparison with Other Codes
- **vs tinydft**: python_1d_dft is even simpler (1D vs 3D atoms) and focuses purely on the algorithm flow.
- **vs PyDFT**: PyDFT handles 3D Gaussian basis; python_1d_dft is real-space 1D.
- **Unique strength**: absolute minimal barrier to entry for understanding the "self-consistent loop".

## Application Areas
- **Classroom Teaching**: Perfect for a single-lecture demo.
- **Self-Study**: For students learning the Kohn-Sham equations.
- **Algorithm Prototyping**: Testing simple functionals or mixing schemes in 1D.

## Best Practices
- **Read the Code**: The source code *is* the documentation.
- **Vary Parameters**: Experiment with electron count and grid density to see effects.
- **Plot Results**: Use the plotting functions to visualize how the density changes during SCF.

## Community and Support
- **GitHub**: Open source repository with issues/discussions.
- **Tutorials**: The repo itself is structured as a tutorial.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/tamuhey/python_1d_dft

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Purpose: Clearly educational and functional for its scope.
