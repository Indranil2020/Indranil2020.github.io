# BERTHA

## Official Resources
- Source Repository: https://github.com/BERTHA-4c-DKS/pybertha
- Documentation: Integrated within repository and PyBERTHA examples.
- License: LGPL-3.0

## Overview
BERTHA is a state-of-the-art relativistic Density Functional Theory (DFT) code. It implements the full four-component Dirac-Kohn-Sham (DKS) formalism, enabling highly accurate electronic structure calculations for systems containing heavy elements where relativistic effects are dominant. The code has recently evolved to include Python bindings (PyBERTHA), making it accessible for modern scripting workflows.

**Scientific domain**: Relativistic quantum chemistry, heavy element chemistry, f-block elements
**Target user community**: Researchers involving actinides, lanthanides, and heavy-element molecular physics

## Theoretical Methods
- **Dirac-Kohn-Sham (DKS) DFT**: Full four-component relativistic treatment.
- **Basis Sets**: G-spinor basis sets (Gauzin-spinors).
- **Hamiltonians**: 
    - Dirac-Coulomb (DC)
    - Dirac-Coulomb-Gaunt (DCG) (interaction treated variationally or perturbatively)
- **Density Fitting**: Efficient evaluation of 4-center integrals using J-fitting.
- **Time-Dependent DFT**: Relativistic TDDFT capabilities for excited states.

## Advanced Features
- **Parity Violation**: 
  - Specialized capabilities for calculating parity-violation energy shifts in chiral molecules.
  - Evaluation of weak nuclear interaction effects.
- **G-Spinor Basis**:
  - Use of kinetically balanced Gaussian-spinor basis sets for optimal relativistic description.
  - Automatic generation of small component basis functions.

## Capabilities
- **Relativistic Electronic Structure**: Accurate ground state energies and properties for heavy atoms and molecules.
- **Python Interface**: `pybertha` allows for flexible input generation, execution, and post-processing.
- **Molecular Properties**: 
    - Hyperfine structure constants
    - Electric field gradients


## Key Strengths
### Rigorous Relativistic Formulation
- Avoids approximations inherent in scalar-relativistic or two-component methods.
- Essential for verifying results of approximate methods on heavy systems.

### Modern Interface
- The PyBERTHA interface integrates the legacy Fortran core with a modern Python environment, facilitating easier workflow automation.

## Inputs & Outputs
- **Input**: Python scripts via `pybertha` or legacy input files.
- **Output**: Energies, orbital coefficients, property values, formatted text output.

## Interfaces & Ecosystem
- **Languages**: Core in Fortran, Interface in Python.
- **Parallelization**: OpenMP for shared-memory parallelism.
- **Dependencies**: Standard scientific python stack (numpy, scipy), Fortran compiler.

## Computational Cost
- **Complexity**: Higher cost than non-relativistic DFT due to 4-component spinor formalism.
- **Optimization**: Uses density fitting and efficient integral algorithms to mitigate the $N^4$ scaling.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/BERTHA-4c-DKS/pybertha
2. "BERTHA: A relativistic density functional theory code" (Search for associated publications)

## Community and Support
- **Support Strategy**: Direct interaction via GitHub Issues.
- **Development**: Active contributions from original authors (Univ. of Chieti-Pescara).

**Confidence**: VERIFIED
**Status**: Active, Open Source
**Note**: Ensure to use the `pybertha` repository for the modern Python interface.
