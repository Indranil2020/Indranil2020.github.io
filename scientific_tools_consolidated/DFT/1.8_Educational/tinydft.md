# tinydft

## Official Resources
- **Homepage**: https://github.com/theochem/tinydft
- **Source Repository**: https://github.com/theochem/tinydft
- **License**: MIT License

## Overview
**tinydft** is a minimalist, educational Density Functional Theory code written in Python. It is designed to be as simple as possible to illustrate the core concepts of DFT without the complexity of production-level optimizations. It focuses on spherically symmetric atoms and serves as a pedagogical tool.

**Scientific domain**: Education, Atomic Physics.
**Target user community**: Students learning DFT, Instructors.

## Theoretical Methods
- **System**: Spherically symmetric atoms (1D radial grid).
- **Functionals**: LDA (Dirac exchange).
- **Discretization**: Finite differences on a radial grid.

## Capabilities
- **Simplicity**: Codebase is extremely small and readable.
- **Concept Demonstration**: Shows how the Kohn-Sham equations are solved iteratively.

## Key Strengths
- **Legibility**: intended to be read like a textbook example.
- **Lightweight**: Minimal dependencies (NumPy).
## Performance Characteristics
- **Scope**: Very fast for small atoms; performance is not a priority compared to code readability.
- **Implementation**: Pure Python/NumPy.

## Limitations & Known Constraints
- **Atoms Only**: Restricted to spherically symmetric atoms; cannot simulate molecules or crystals.
- **Functionals**: Default support is limited to Dirac exchange (LDA).
- **Physics**: Uses rigid occupation rules (Klechkowski) and spin-unpolarized densities by default (though extensible).

## Best Practices
- **Usage**: Use the `tinydft-demo-mendelejev` script to visualize radial densities immediately.
- **Extension**: The best way to use the code is to attempt the "programming assignments" provided in the repo (e.g., adding Hartree-Fock).

## Community and Support
- **Support**: Educational tool provided "as-is".
- **Resources**: The repository acts as a self-contained tutorial.
## Verification & Sources
**Primary sources**:
1.  **Repository**: [tinydft GitHub](https://github.com/theochem/tinydft)

**Verification status**: ✅ VERIFIED
