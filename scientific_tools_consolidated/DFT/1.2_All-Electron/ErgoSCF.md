# ErgoSCF

## Official Resources
- Homepage: http://www.ergoscf.org/
- Source Repository: https://github.com/ergoscf/ergoscf  (or SourceForge)
- Documentation: http://www.ergoscf.org/documentation.html
- License: GNU General Public License (GPL)

## Overview
ErgoSCF is a quantum chemistry program designed for large-scale, linear-scaling electronic structure calculations. It works with Gaussian basis sets and enforces a strict all-electron methodology. The code is built for efficiency, utilizing modern techniques like fast multipole methods and sparse matrix algebra to handle large molecules and clusters.

**Scientific domain**: Large molecules, clusters, rigorous all-electron chemistry
**Target user community**: Quantum chemists needing linear-scaling all-electron calculations

## Theoretical Methods
- Density Functional Theory (DFT) (Kohn-Sham)
- Hartree-Fock (HF)
- Linear Scaling (O(N)) methodology
- Gaussian basis sets
- All-electron formulations (no effective core potentials by default)
- Fast Multipole Method (FMM) for Coulomb interactions
- Sparse matrix algebra

## Capabilities
- Ground-state energy and gradients
- Geometry optimization
- Linear response properties (polarizabilities)
- Calculation of large molecular systems
- Harmonic vibrational frequencies
- Population analysis

## Key Strengths
### Linear Scaling
- Achieves O(N) scaling for both HF and DFT
- Enables all-electron calculations on large systems

### Efficiency
- Hierarchical sparse matrix infrastructure
- Rigorous integral screening
- Memory efficient for large basis sets

## Inputs & Outputs
- **Input**: Text-based input files describing geometry and calculation parameters
- **Output**: Standard output with energies, properties, and optimized geometries

## Interfaces & Ecosystem
- **Programming Language**: C++
- **Python Integration**:
  - No native Python API.
  - Can be orchestrated via `mpi4py` or `multiprocessing` for workflows.
- **Parallelization**:
  - Hybrid MPI/OpenMP support for distributed and shared memory systems.

## Advanced Features
- **Hierarchical Matrices**:
  - Exploits sparsity in large systems for efficiency.
- **Fast Multipole Method (FMM)**:
  - Accelerates Coulomb interaction calculations to O(N).
- **Rigorous Screening**:
  - Exact control over error bounds, ensuring reliability.

## Community and Support
- **Documentation**: Available at http://www.ergoscf.org/documentation.html
- **Support**: Via website contacts and repository issues.

## Computational Cost
- **Scaling**: Linearly scaling O(N) with system size for large systems.
- **Memory**: O(N) memory usage, allowing calculations on commodity hardware.

## Verification & Sources
**Primary sources**:
1. Official Website: http://www.ergoscf.org/
2. "Ergo: An open-source program for linear-scaling electronic structure calculations"

**Confidence**: VERIFIED
**Status**: Mature, Open Source
