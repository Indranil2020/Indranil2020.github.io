# QMC2 (jorgehog)

## Official Resources
- Homepage: https://github.com/jorgehog/QMC2
- Source Repository: https://github.com/jorgehog/QMC2
- License: Open Source (Implicit/Check Repo)

## Overview
**QMC2** is a C++ implementation of quantum Monte Carlo methods, specifically focusing on an efficient **Diffusion Monte Carlo (DMC)** algorithm. Developed by *jorgehog*, it combines a high-performance C++ core with flexible Python scripting for setup and analysis. It is designed for calculating ground state properties of atoms, molecules, and extended systems with high precision.

**Scientific domain**: Quantum Chemistry, Electronic Structure, Many-Body Physics
**Target user community**: Developers and researchers in QMC methods

## Theoretical Methods
- **Diffusion Monte Carlo (DMC)**: Real-space fixed-node DMC for projecting ground states.
- **Variational Monte Carlo (VMC)**: Optimization of trial wavefunctions.
- **Trial Wavefunctions**: Supports Slater-Jastrow type wavefunctions.
- **Basis Sets**: Implements analytical potentials and basis functions (Gaussian/Slater).

## Capabilities (CRITICAL)
- **Efficient Solver**: Optimized C++ implementation of the DMC random walk.
- **Parallelization**: MPI-based parallelism for scaling across multiple cores/nodes.
- **Custom Potentials**: Extensible class structure (`Potentials` subclassing) to define new physical systems.
- **Orbital Management**: `OrbitalGenerator` (via SymPy) and `Orbitals` classes for handling complex basis sets and mappings.
- **Visualization**: Integration points for `DCViz` and PySide for observing walkers/densities (optional GUI components).

## Key Features

### Hybrid Architecture:
- **C++ Core**: Heavy numerical lifting (random walks, wavefunction evaluation) in compiled C++.
- **Python Layer**: Simulation setup, parameter definition, and result analysis in Python.

### External Libraries:
- **Armadillo**: For efficient linear algebra operations.
- **Boost**: For general utilities and mathematical functions.
- **MPI**: For parallel implementation (Optional).

### Extensibility:
- designed for users to plug in new hamiltonians or basis sets without rewriting the core engine.

## Inputs & Outputs
- **Input formats**:
  - Python scripts defining the system geometry and trial function parameters.
  - XML/JSON configurations (if applicable).
- **Output data types**:
  - Text/Binary logs of energies, variance, and acceptance ratios.
  - Checkpoint files for restarting walkers.

## Interfaces & Ecosystem
- **Dependencies**: Armadillo, Boost, MPI, Python headers.
- **Visualization**: Connects with `DCViz` for data visualization.

## Workflow and Usage
Users typically compile the C++ backend and then write a Python script that imports the QMC2 module. The script defines the atoms, basis set, and QMC parameters (time step, walkers), then invokes the VMC optimization followed by the DMC propagation.

## Performance Characteristics
- **Speed**: Efficient C++ memory management and math kernels.
- **Scaling**: Linearly scalable with number of MPI processes (walker parallelism).

## Comparison with Other Codes
| Feature | QMC2 | QMCPACK | CASINO |
| :--- | :--- | :--- | :--- |
| **Scope** | Research / Specialized DMC | General Production / HPC | General Production / Molecules |
| **Language** | C++ Core / Python Scripting | C++ / CUDA | Fortran 95 |
| **Performance** | Good (MPI) | Excellent (Exascale/GPU) | High |
| **Ecosystem** | Smaller / Custom | Large / Community Standard | Large / Academic |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/jorgehog/QMC2

**Verification status**: ✅ VERIFIED
- Source code: OPEN
- Status: Research code, active in past years.
- Focus: Efficient C++ DMC implementation.
