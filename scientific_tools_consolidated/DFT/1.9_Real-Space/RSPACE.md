# RSPACE

## Official Resources
- Homepage: https://materiapps.issp.u-tokyo.ac.jp/en/apps/rspace/
- Distribution: Part of MateriApps suite
- License: Contact developers / Proprietary/Academic functionality

## Overview
RSPACE is a first-principles simulation code package based on the real-space finite-difference method using pseudopotentials. It is specifically tailored for high-speed and high-precision calculations of electronic states in aperiodic systems such as surfaces, solid interfaces, clusters, and nanostructures. It offers specialized implementation for quantum transport properties under semi-infinite boundary conditions.

**Scientific domain**: Surface science, Quantum transport, Nanostructures
**Target user community**: Researchers in spintronics, transport, and surface physics

## Theoretical Methods
- Real-Space Finite-Difference Method (RSFD)
- Density Functional Theory (DFT)
- Projector Augmented Wave (PAW) method
- Norm-Conserving Pseudopotentials
- Overbridging Boundary Matching (OBM) method for transport
- Non-equilibrium Green's Function (NEGF) - related methods

## Capabilities
- **Electronic Structure**: Band structures, DOS for large supercells.
- **Quantum Transport**: Conductance calculations for intrinsic nanostructures.
- **Magnetism**: Spin-orbit interaction, non-collinear magnetism.
- **Boundary Conditions**: 0D (cluster), 1D (wire), 2D (film/surface), 3D (bulk).

## Key Strengths

### Transport Calculations:
- Specialized for calculating electron transport through nanostructures bridging semi-infinite electrodes.
- Efficient handling of open boundary conditions.

### Versatile Grid Method:
- No basis set superposition error (BSSE).
- Flexible boundary conditions without vacuum padding issues.

### PAW Implementation:
- Accurate treatment of transition metals and magnetic systems using PAW in real-space.

## Inputs & Outputs
- **Inputs**:
  - Grid parameters
  - Structure file
  - Transport boundary definitions
- **Outputs**:
  - Transmission coefficients
  - Current-voltage (I-V) characteristics
  - Spin-resolved densities

## Interfaces & Ecosystem
- **MateriApps**: Integrated into the MateriApps Live! environment.
- **Visualization**: Output compatible with VESTA and other standard tools via conversion.

## Advanced Features
- **Krylov Subspace**: Efficient iterative solvers for large sparse matrices.
- **Spin Dynamics**: Non-collinear spin texture analysis.

## Performance Characteristics
- **Parallelization**: Parallelized via MPI domain decomposition.
- **Scalability**: High scalability due to locality of finite difference operators.

## Computational Cost
- **Moderate to High**: Transport calculations are computationally demanding; real-space grid requires fine meshing for deep potentials.

## Limitations & Known Constraints
- **Availability**: Distribution seems less "open" than standard GitHub repos; often obtained via MateriApps or direct contact.
- **Documentation**: Primary resources are often in Japanese or technical reports; English documentation varies.

## Comparison with Other Codes
- **vs TranSIESTA**: TranSIESTA uses LCAO basis; RSPACE uses real-space grid (more accurate but more costly).
- **vs OpenMX**: OpenMX (LCAO) is also strong in transport; RSPACE offers a basis-set-free alternative check.
- **Unique strength**: Real-space formulation of quantum transport with PAW accuracy.

## Application Areas
- **Spintronics**: Magnetic tunnel junctions, spin filters.
- **Molecular Electronics**: Single-molecule junctions.
- **Surface Reactions**: Catalysis on surfaces (aperiodic).

## Best Practices
- **Grid Sizing**: Ensure grid is fine enough for PAW projectors.
- **Transport**: carefully define electrode regions.

## Community and Support
- **Origin**: Developed by groups at Osaka University, University of Tsukuba, and others.
- **Support**: via MateriApps forums.

## Verification & Sources
**Primary sources**:
1. MateriApps Profile: https://materiapps.issp.u-tokyo.ac.jp/en/apps/rspace/
2. K. Hirose et al., "First-Principles Calculations in Real-Space Formalism", Imperial College Press (2005).

**Verification status**: ✅ VERIFIED
- Existence: Confirmed via MateriApps and publications.
- Accessibility: Available via specific academic channels.
