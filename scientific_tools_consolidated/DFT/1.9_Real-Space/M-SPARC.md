# M-SPARC

## Official Resources
- Homepage: https://github.com/SPARC-X/M-SPARC
- Source Repository: https://github.com/SPARC-X/M-SPARC
- License: Open Source (GPLv3)

## Overview
M-SPARC is the MATLAB version of the SPARC (Simulation Package for Ab-initio Real-space Calculations) electronic structure code. It provides a rapid prototyping environment for real-space density functional theory (DFT) calculations. Leveraging MATLAB's powerful linear algebra capabilities, M-SPARC offers an accessible codebase for developing new algorithms and understanding the real-space finite-difference method without the complexity of C++/MPI implementations.

**Scientific domain**: Real-Space DFT, Prototyping, Education
**Target user community**: Algorithm developers, Students, Researchers prototyping new functional/methods

## Theoretical Methods
- Real-Space Finite-Difference Method
- Kohn-Sham Density Functional Theory
- Pseudopotentials (Norm-conserving)
- Chebyschev Filtering (for eigensolver)
- Newton-Raphson Optimization
- Periodic and Non-Periodic Boundary Conditions

## Capabilities
- **Electronic Structure**: Ground state energy, forces, band structure.
- **Molecular Dynamics**: Ab initio MD (NVE, NVT).
- **Relaxation**: Geometry optimization of atomic positions.
- **Boundary Conditions**: Isolated, Periodic, Mixed.
- **Spin**: Spin-polarized and unpolarized calculations.
- **Symmetry**: Cyclical symmetry exploitation.

## Key Strengths

### Rapid Prototyping:
- MATLAB implementation allows for extremely fast implementation of new ideas.
- Code structure mirrors the mathematical formulation of Real-Space DFT.

### Real-Space Precision:
- High-order finite difference stencils for systematic convergence.
- No FFTs required; purely local or sparse-matrix operations.

### Accessibility:
- Significantly easier to read and modify than the main C++ SPARC code.
- Ideal for teaching the details of real-space DFT implementation.

## Inputs & Outputs
- **Inputs**:
  - MATLAB script/structure defining parameters (system, grid, mixing).
  - Pseudopotential files (psp8, upf).
  - Atomic coordinates.
- **Outputs**:
  - Energy components.
  - Forces.
  - Charge density plots (MATLAB figures or cube/xsf export).
  - Band structures.

## Interfaces & Ecosystem
- **MATLAB**: Full integration with MATLAB toolbox.
- **SPARC**: Serves as the development/testing ground for features eventually ported to the C++ SPARC code.

## Advanced Features
- **Linear Scaling**: Implementation of Chebyschev Filtering methods.
- **Krylov Solvers**: Access to MATLAB's robust iterative solvers.

## Performance Characteristics
- **Speed**: Slower than optimal C++ codes due to MATLAB overhead, but efficient for small-to-medium systems.
- **Parallelism**: Uses MATLAB's native multi-threading.

## Computational Cost
- **Moderate**: Best for unit cells up to ~100 atoms or prototyping; for large-scale production, use C++ SPARC.

## Limitations & Known Constraints
- **License Requirement**: Requires proprietary MATLAB license.
- **Performance**: Not suitable for massive HPC runs like the C++ version.

## Comparison with Other Codes
- **vs SPARC (C++)**: M-SPARC is the development/teaching twin; SPARC is the production HPC twin.
- **vs KSSOLV**: Both are MATLAB DFT codes; M-SPARC focuses on Real-Space Finite-Difference, KSSOLV on Plane-Wave/Real-Space mix.
- **Unique strength**: The cleanest entry point for learning Real-Space DFT algorithms.

## Application Areas
- **Algorithm Development**: Testing new solvers, mixing schemes, or functionals.
- **2D Materials**: Testing boundary conditions for slabs/sheets.
- **High Temperature**: Warm dense matter method testing.

## Best Practices
- **Use for Dev**: Write new features here, debug them, then port to C++.
- **Grid Check**: Always verify `MESH_SPACING` convergence.
- **Pseudopotentials**: Ensure compatibility with SPARC-formatted potentials.

## Community and Support
- **Source**: Maintained by SPARC-X group (Georgia Tech, etc.).
- **Response**: GitHub issues active.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/SPARC-X/M-SPARC
2. Publication: "M-SPARC: A MATLAB implementation of the simulation package..." SoftwareX (2020).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GPLv3) but requires MATLAB.
