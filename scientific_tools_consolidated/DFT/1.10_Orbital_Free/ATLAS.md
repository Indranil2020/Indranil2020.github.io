# ATLAS (Orbital-Free DFT)

## Official Resources
- Repository/Source: Research code (referenced in literature, specific repo private/restricted)
- Key Publication: *Mi et al., Comp. Phys. Comm. 200, 87 (2016)*
- License: Academic/Research

## Overview
ATLAS is a specialized Orbital-Free Density Functional Theory (OF-DFT) software package. Unlike Kohn-Sham DFT, ATLAS solves for the electron density directly without using orbitals, typically employing a real-space finite-difference method. This allows it to scale linearly with system size (O(N)) with a very small prefactor, enabling simulations of millions of atoms, particularly for main-group metals (like Al, Mg) where OF-DFT kinetic energy functionals are accurate.

**Scientific domain**: Large-scale metallurgy, Warm Consensus Matter, High-pressure physics
**Target user community**: Method developers, Metallurgists modeling dislocations/grain boundaries

## Theoretical Methods
- Orbital-Free Density Functional Theory (OF-DFT)
- Real-Space Finite-Difference method
- Wang-Teter / Thomas-Fermi-von Weizsäcker Kinetic Energy Functionals
- Local pseudopotentials (bulk-derived)
- Non-local density kernels (for kinetic energy)
- Newton-Krylov optimization algorithms

## Capabilities
- **System Size**: Routine >100,000 atoms; Capable of millions.
- **Materials**: Primarily simple metals (Al, Mg, Li, Na) and their alloys.
- **Output**: Direct equilibrium density, total energy, forces, stress.
- **Simulation**: Geometry optimization, Molecular Dynamics.

## Key Strengths

### Linear Scaling O(N):
- By avoiding orbitals and orthogonalization, computational cost scales purely linearly with volume/atoms.
- Allows simulation of mesoscopic features (dislocations, voids) at quantum accuracy.

### Robust Minimization:
- Implements advanced numerical algorithms (trust-region Newton methods) to ensure convergence of the highly non-linear OF-DFT energy functional.

## Inputs & Outputs
- **Inputs**:
  - Cell dimensions
  - Ion positions
  - Functional choice (TF, vW, WT, etc.)
  - Local Pseudopotential parameters
- **Outputs**:
  - Total Energy
  - Electron Density (real space grid)
  - Forces

## Interfaces & Ecosystem
- **Standalone**: Typically runs as a standalone solver.
- **Integration**: Methods often cross-pollinate with PROFESS.

## Advanced Features
- **Mixed Basis**: Some versions explore mixed basis or adaptive grids.
- **Parallelization**: MPI parallelization over spatial domains.

## Performance Characteristics
- **Speed**: Extremely fast compared to KS-DFT for large systems (orders of magnitude).
- **Memory**: very low memory footprint (density only, no wavefunctions).

## Computational Cost
- **Very Low**: Can simulate thousands of atoms on a workstation.

## Limitations & Known Constraints
- **Accuracy**: Limited by the accuracy of the Kinetic Energy Density Functional (KEDF). High accuracy only for "simple" metals (Al, Mg). Failed for Transition Metals (d-electrons) and covalent bonds (Si, C) without advanced/modern functionals.
- **Pseudopotentials**: Must use Local Pseudopotentials (LPS), which transfers less well than Non-Local ones.

## Comparison with Other Codes
- **vs PROFESS**: PROFESS is the leading open-source OF-DFT code; ATLAS is a robust alternative often used for specific benchmarks or internal research.
- **vs VASP**: ATLAS is fundamentally different (orbital-free vs Kohn-Sham); ATLAS wins on size, VASP wins on versatility/accuracy.
- **Unique strength**: Specialized optimization algorithms for solving the Euler-Lagrange equation of OF-DFT.

## Application Areas
- **Metallurgy**: Study of grain boundaries, dislocations, and stacking faults in Aluminum/Magnesium.
- **Liquid Metals**: Large scale MD of liquid metals.
- **Warm Dense Matter**: High temperature electron gas simulations.

## Best Practices
- **Check Material**: Only apply to Al/Mg/Li unless you are testing new functionals.
- **Validate KEDF**: Ensure the chosen Kinetic Energy Functional reproduces KS-DFT properties for the bulk phase first.

## Community and Support
- **Research Based**: Support is via the authors of the key papers (e.g., W. Mi, H. Wang).

## Verification & Sources
**Primary sources**:
1. Mi, W., et al., "ATLAS: A real-space finite-difference implementation of orbital-free density functional theory", *Computer Physics Communications*, 200, 87-95 (2016).

**Verification status**: ✅ VERIFIED
- Existence: Confirmed via widely cited CPC publication.
- Availability: Research-grade distribution.
