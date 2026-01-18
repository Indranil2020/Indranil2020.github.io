# PyTreeNet

## Official Resources
- Homepage: https://github.com/Drachier/PyTreeNet
- Documentation: arXiv:2407.13249 ("PyTreeNet: A Python Library for easy Utilisation of Tree Tensor Networks")
- Source Repository: https://github.com/Drachier/PyTreeNet
- License: EUPL v1.2

## Overview
**PyTreeNet** is a Python library dedicated to the simulation of quantum many-body systems using **Tree Tensor Networks (TTN)**. Developed by the "Drachier" team, it generalizes Matrix Product States (MPS) to tree-like topologies, allowing for efficient representation of systems with hierarchical entanglement structures or non-1D controlivities. It focuses on easing the implementation of complex tensor network algorithms like ground state search and time evolution.

**Scientific domain**: Tensor Networks, Quantum Information, Many-Body Physics
**Target user community**: Researchers experimenting with TTN algorithms and hierarchical quantum systems

## Theoretical Methods
- **Tree Tensor Networks (TTN)**: Acyclic graph states generalizing MPS.
- **Time Evolution**:
  - **TEBD**: Time-Evolving Block Decimation generalized to trees.
  - **TDVP**: Time-Dependent Variational Principle for optimal manifold projection.
- **Ground State Search**: Variational algorithms akin to DMRG but on tree structures.

## Capabilities (CRITICAL)
- **Arbitrary Trees**: Supports defining Hamiltonian and States on arbitrary tree graphs.
- **Symbolic Input**: Converts symbolic Hamiltonians into TTN operator forms (Tree MPO).
- **Dynamics**: Native support for real-time evolution, crucial for studying non-equilibrium quench dynamics.
- **Pythonic API**: high-level abstraction built on NumPy, hiding the complexity of tensor index contractions.

## Key Features

### Usability:
- **Easy Installation**: `pip install pytreenet`.
- **Documentation**: Accompanied by a detailed arXiv guide with exercises and examples.

### Algorithms:
- **TDVP Integration**: Implements the robust TDVP algorithm, often superior to TEBD for long-time evolution and Hamiltonians with long-range terms.

## Inputs & Outputs
- **Input formats**:
  - Python scripts defining the graph connectivity and local Hilbert spaces.
  - Symbolic Hamiltonian expressions.
- **Output data types**:
  - Expectation values of observables.
  - Evolved states.

## Interfaces & Ecosystem
- **Dependencies**: NumPy, SciPy.
- **Ecosystem**: Can be used alongside other Python TN libraries, but functions as a standalone high-level tool.

## Workflow and Usage
```python
import pytreenet as ptn
# Define tree structure and Hamiltonian
tree = ptn.Tree(...)
H = ptn.Hamiltonian(...)
# Ground state search
psi = ptn.dmrg(H, tree)
# Time evolution
psi_t = ptn.tdvp(psi, H, dt=0.01, steps=100)
```

## Performance Characteristics
- **Efficiency**: Good for systems where entanglement follows a tree hierarchy (e.g., dendrimers, effective models).
- **Backend**: NumPy-based (CPU), suitable for prototyping and intermediate scale problems.

## Comparison with Other Codes
| Feature | PyTreeNet | ITensor | TeNPy |
| :--- | :--- | :--- | :--- |
| **Network Type** | Tree Tensor Network (TTN) | MPS (1D) / MPO | MPS (1D) |
| **Language** | Python (NumPy) | C++ / Julia | Python |
| **Connectivity** | Hierarchical / Arbitrary Tree | Linear (mostly) | Linear |
| **Focus** | Complex topologies / Dendrimers | 1D Physics | 1D Physics |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/Drachier/PyTreeNet
2. Publication: "PyTreeNet: A Python Library for easy Utilisation of Tree Tensor Networks", arXiv:2407.13249.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (EUPL v1.2)
- State: Active research code, documented in 2024.
