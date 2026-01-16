# S/PHI/nX

## Official Resources
- Homepage: https://sxlib.mpie.de/ (or https://www.mpie.de/2586717/sphinX)
- Documentation: https://sxlib.mpie.de/documentation.html
- Source Repository: https://github.com/SPHInX-X/sxaccelerate (Core library) / Binaries via Open Build Service
- License: Apache License 2.0 / LGPL (varies by module)

## Overview
S/PHI/nX is a C++ based library and software package for electronic structure theory, developed at the Max-Planck-Institut für Eisenforschung. It combines standard plane-wave pseudopotential density functional theory (DFT) with k.p theory and other specialized methods. It is built upon the SxAccelerate library, emphasizing modularity, efficient memory handling, and modern C++ design.

**Scientific domain**: Materials science, thermodynamics, defect physics, surfaces
**Target user community**: Researchers in computational materials science, particularly those needing specialized defect treatments or k.p methods

## Theoretical Methods
- Density Functional Theory (DFT)
- Planewave basis sets
- Norm-conserving and PAW (Projector Augmented Wave) potentials
- k.p perturbation theory
- Hubbard U corrections (DFT+U)
- Hybrid functionals
- Van der Waals corrections
- Ab initio thermodynamics

## Capabilities
- Ground-state electronic structure
- Geometry optimization (BFGS, substitutions)
- Molecular Dynamics (MD)
- Charged defect calculations (sxdefectalign)
- Band structure and Density of States (DOS)
- Stress and force calculations
- Surface and 2D material simulations
- Optical properties

## Key Strengths

### Modern C++ Architecture:
- Built on SxAccelerate for high-performance I/O and data management
- Modular and extensible design
- Object-oriented structure

### Defect Physics:
- Specialized tools for charged defects (sxdefectalign)
- Corrections for finite-size errors in supercells
- 2D material defect alignment (sxdefectalign2d)

### Integration:
- Integrated with pyiron (Python-based IDE for materials science)
- Flexible hierarchical input format (sx format)

## Inputs & Outputs
- **Input formats**:
- **Input formats**:
  - `input.sx`: Hierarchical, block-structured input/output format (similar to JSON/C-structs).
  - Example: `structure { species { ... } }`
  - Associative array style allows flexible parameter definition.
  - Pseudopotentials (standard formats)
  - Structure files
  
- **Output data types**:
  - Energy, Forces, Stress
  - Wavefunctions
  - Charge densities
  - Band structures
  - Thermodynamic data

## Interfaces & Ecosystem
- **Python**: Native integration with **pyiron** for high-level workflow management.
- **SxAccelerate**: Core library available for custom tool development.
- **Visualization**: Output compatible with standard visualization tools (e.g., VESTA via conversion).

## Performance Characteristics
- **Speed**: C++ optimized performance.
- **Parallelizaton**: MPI parallelization for large-scale runs.
- **Efficiency**: Efficient memory handling via SxAccelerate pointers.

## Best Practices

### Input File Management:
- **Modularity**: Use the hierarchical format to organize complex simulations.
- **Comments**: Heavily comment `input.sx` files (C++ style `//`) for reproducibility.

### Defect Calculations:
- **Charged Defects**: Always use `sxdefectalign` to correct for finite-size electrostatic errors in supercells.
- **Relaxation**: Use BFGS for robust geometry optimization of defect structures.

### Workflow:
- **pyiron**: Leverage the `pyiron` IDE to manage job submission and database storage of results.
- **Compilation**: Use the provided binaries or containers if compiling from source proves difficult due to C++ dependencies.

## Community and Support
- **Hosting**: Active on [GitHub](https://github.com/SPHInX-X).
- **Organization**: Maintained by Max-Planck-Institut für Eisenforschung (MPIE).
- **Documentation**: Comprehensive C++ API docs and user tutorials online.
- **Issues**: Issue tracking via GitHub repository.

## Limitations & Known Constraints
- **Community**: Smaller user base compared to VASP or Quantum ESPRESSO.
- **Complexity**: The flexible input format and C++ structure generally have a learning curve.
- **Compilation**: S/PHI/nX compilation can be complex due to dependencies, though binaries are provided.

## Comparison with Other Codes
- **vs VASP**: Similar plane-wave capabilities; S/PHI/nX offers open-source freedom and specialized C++ structure, but VASP has a broader feature set.
- **vs Quantum ESPRESSO**: Both are open source; S/PHI/nX focuses more on unique C++ modularity and defect physics features.

## Verification & Sources
**Primary sources**:
1. Max-Planck-Institut website: https://www.mpie.de/2586717/sphinX
2. SxAccelerate GitHub: https://github.com/SPHInX-X
3. Freysoldt et al., "Fully self-consistent GW calculations..." (S/PHI/nX usage in literature)

**Confidence**: CONFIRMED - Established code from a reputable institute (MPIE).

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: DFT/Plane-Wave
- Key Feature: C++ Library, Defect Physics
