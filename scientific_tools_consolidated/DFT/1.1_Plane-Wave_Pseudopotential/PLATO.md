# PLATO (Package for Linear Combination of Atomic Orbitals)

## Official Resources
- **Homepage (Legacy)**: http://www.dl.ac.uk/TCSC/Software/PLATO/
- **Research Group**: https://www.imperial.ac.uk/materials/research/tsm/ (Theory and Simulation of Materials, Imperial College London)
- **References**: S. D. Kenny, A. P. Horsfield, H. Fujitani, Phys. Rev. B 62, 4899 (2000).
- **License**: Academic/Research (typically distributed via CCP9 or request)

## Overview
PLATO is a localized orbital-based electronic structure package developed by Andrew Horsfield, Steve Kenny, and collaborators (Imperial College London, UCL, Loughborough). It allows for both tight-binding and density functional theory (DFT) calculations within a single framework. PLATO is particularly noted for its efficiency in handling large systems using O(N) methods and its versatility in treating both orthogonal and non-orthogonal basis sets.

**Scientific domain**: Tight-binding, DFT, localized orbitals, O(N) methods  
**Target user community**: Materials scientists, tight-binding researchers, large-scale simulation community

## Theoretical Methods
- Density Functional Theory (DFT)
- Tight-Binding (TB)
- Numerical Atomic Orbitals (NAO)
- Linear Combination of Atomic Orbitals (LCAO)
- Non-orthogonal and orthogonal basis sets
- Multipole expansions for electrostatics
- O(N) scaling algorithms

## Capabilities
- Ground-state electronic structure
- Structural relaxation
- Molecular dynamics
- Tight-binding parametrization
- DFT calculations with localized bases
- Large system simulations (thousands of atoms)
- Transport properties (with additional modules)
- Point defects and extended defects

## Key Strengths
### Unified Framework
- Seamlessly bridges Tight-Binding and DFT
- Allows testing of TB parameters against DFT within the same code

### Efficiency
- Optimized for O(N) scaling
- Efficient handling of large supercells
- Suitable for complex defect structures

## Inputs & Outputs
- **Input formats**:
  - `input.dat` (Control parameters)
  - Structure files
  - Basis set definitions
- **Output data types**:
  - Energy, forces, stress
  - Charge densities
  - Band structures (using auxiliary tools)

## Limitations
- **Availability**: Not a standard open-source repository; often obtained via academic channels (CCP9).
- **Documentation**: Less publicly accessible than major community codes like VASP or QE.

## Computational Cost
- **O(N) Scaling**: Linear scaling with system size, enabling calculations on thousands of atoms.
- **Efficiency**: Very high for tight-binding; DFT mode slower but competitive for large systems.

## Comparison with Other Codes
- **vs SIESTA**: Both are O(N) codes using localized orbitals. PLATO allows direct comparison of TB and DFT parameters.
- **vs VASP**: PLATO is specialized for large systems/O(N); VASP is a general purpose plane-wave code (scaling $N^3$).

## Best Practices
- **Basis Set**: Careful testing of basis set completeness is required (unlike plane waves).
- **Parameters**: TB parameters must be transferable to your system of interest.

## Community and Support
- **Access**: Distributed via CCP9 (Collaborative Computational Project 9) or by request.
- **Support**: Direct academic collaboration with developers.

## Verification & Sources
- **Primary Source**: Published literature (Phys. Rev. B 62, 4899) and Imperial College research pages.
- **Confidence**: VERIFIED - Well-established code in the tight-binding/DFT community.
