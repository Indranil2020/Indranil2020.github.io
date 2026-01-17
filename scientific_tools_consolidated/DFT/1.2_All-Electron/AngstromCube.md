# AngstromCube

## Official Resources
- Source Repository: https://github.com/real-space/AngstromCube
- Documentation: https://github.com/real-space/AngstromCube (README)
- License: Open Source (Check repository for specific license)

## Overview
AngstromCube is a high-performance Density Functional Theory (DFT) code designed for large-scale all-electron calculations. It employs a real-space finite-difference formulation, allowing for efficient parallelization and O(N) linear scaling. The code is specifically optimized for modern hardware, featuring native GPU acceleration to handle systems with thousands of atoms.

**Scientific domain**: Materials science, large-scale nanostructures, all-electron precision
**Target user community**: Researchers needing high-accuracy all-electron calculations for large systems, HPC users

## Theoretical Methods
- Density Functional Theory (DFT)
- Real-Space Finite Difference method
- All-Electron methodology (no pseudopotentials)
- Linear-scaling O(N) algorithms
- Real-time Time-Dependent DFT (RT-TDDFT) capability
- Chebyshev filtering for eigensolvers

## Capabilities
- Ground-state electronic structure
- Large-scale system simulations (thousands of atoms)
- Massively parallel execution (MPI + GPU)
- All-electron precision for core states
- Efficient scale-up on supercomputers

## Key Strengths
### Real-Space Grid
- Systematic convergence via grid spacing
- No basis set superposition error (BSSE)
- Flexible boundary conditions

### Performance
- GPU-accelerated (CUDA)
- Linear scaling for large systems
- High parallel efficiency

## Inputs & Outputs
- **Input**: Configuration files (typically text/script based), atomic coordinates
- **Output**: Electronic density, total energy, forces, wavefunctions

## Interfaces & Ecosystem
- **Python Integration**:
  - Native Python interface available for scripting and workflow control.
  - Can be integrated with ASE via custom calculators (not native yet).
- **Data formats**:
  - Supports standard cube file outputs for density/potentials.
  - Compatible with visualization tools like VESTA and VMD via standard formats.

## Advanced Features
- **Real-Time TDDFT**:
  - Simulation of electron dynamics in real-time.
  - Suitable for strong field physics and optical response.
- **Chebyshev Filtering**:
  - Efficient eigensolver for interior eigenvalues.
  - Accelerates convergence for large basis sets.

## Performance Characteristics
- **Speed**: High throughput on GPU systems.
- **Scaling**: O(N) linear scaling with system size.
- **Parallelization**: Hybrid MPI + GPU offloading.

## Community and Support
- **Development**: Active on GitHub.
- **Issues**: Bug tracking via GitHub Issues.

## Computational Cost
- **Scaling**: O(N) for large systems, reducing the cubic scaling bottleneck of traditional DFT.
- **Hardware**: Optimized for GPU architectures, reducing time-to-solution.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/real-space/AngstromCube
2. "AngstromCube: A parallel and GPU-accelerated code for Real-Space All-Electron Linear-Scaling Density Functional Theory"

**Confidence**: VERIFIED
**Status**: Active Development
**Note**: "AngstromCube" also refers to the cubic angstrom unit; ensure searches target the software repository.
