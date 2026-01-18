# DMFTpack

## Official Resources
- Homepage: https://dmftpack.github.io/
- Source Repository: https://github.com/dmftpack/dmftpack
- License: LGPL-3.0

## Overview
DMFTpack is a software package designed for performing DFT+DMFT (Density Functional Theory + Dynamical Mean-Field Theory) calculations. It provides a bridge between first-principles calculations and many-body techniques. A key feature is its inclusion of native impurity solvers such as Iterative Perturbation Theory (IPT) and Self-Consistent Second-Order Perturbation Theory (SC2PT), as well as interfaces to external solvers.

**Scientific domain**: Strongly correlated electron systems, Materials science, Electronic structure
**Target user community**: Researchers in condensed matter physics studying correlated materials

## Theoretical Methods
- Dynamical Mean-Field Theory (DMFT)
- Density Functional Theory (DFT) integration
- Iterative Perturbation Theory (IPT)
- Self-Consistent Second-Order Perturbation Theory (SC2PT)
- Continuous-Time Quantum Monte Carlo (CT-QMC) interface

## Capabilities
- DFT+DMFT self-consistent loop
- Native impurity solvers (IPT, SC2PT)
- Support for external solvers (e.g., CT-QMC)
- Interface with OpenMX for DFT input
- Calculation of spectral functions
- Handling of local interactions

## Key Strengths
### Integrated Solvers:
- Comes with built-in IPT and SC2PT solvers, reducing the need for external dependencies for perturbative approaches.
### OpenMX Interface:
- Specifically designed to interface with the OpenMX DFT code, facilitating efficient workflows.
### Accessibility:
- targeted at making DMFT calculations more accessible with standard perturbative methods.

## Inputs & Outputs
- **Input formats**:
  - OpenMX output files
  - Parameter files for DMFT loop and solver settings
- **Output data types**:
  - Self-energies
  - Green's functions
  - Spectral functions
  - Occupations

## Interfaces & Ecosystem
- **DFT Codes**:
  - OpenMX
- **External Solvers**:
  - CT-QMC (via interface)

## Advanced Features
- **Projectors**: Implements projection methods for mapping DFT Bloch states to local orbitals.
- **Perturbative Approaches**: Efficient solvers for regimes where perturbation theory is valid.

## Performance Characteristics
- **Efficiency**: IPT and SC2PT are computationally much cheaper than QMC methods, allowing for faster calculations on applicable systems.
- **Scalability**: Depends on the underlying impurity solver and the specific system size.

## Computational Cost
- **Low**: For IPT/SC2PT solvers.
- **High**: If using external CT-QMC solvers for complex impurities.

## Limitations & Known Constraints
- **Perturbative Solvers**: IPT and SC2PT are limited to specific parameter regimes (e.g., weak to intermediate coupling) and may fail in strong coupling or specifically for certain multi-orbital problems compared to exact methods like QMC.
- **Documentation**: May be less extensive than major community codes like TRIQS.

## Comparison with Other Codes
- **vs TRIQS**: DMFTpack is a more specific package, whereas TRIQS is a general library for building solvers.
- **vs DFT+DMFT in VASP/QE**: Standalone package focusing on the OpenMX ecosystem initially.

## Verification & Sources
**Primary sources**:
1. Homepage: https://dmftpack.github.io/
2. GitHub: https://github.com/dmftpack/dmftpack

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN
