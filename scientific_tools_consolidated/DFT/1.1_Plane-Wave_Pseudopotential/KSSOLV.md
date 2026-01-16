# KSSOLV

## Official Resources
- Homepage: http://kssolv.org/
- Documentation: http://kssolv.org/ (and included manuals)
- Source Repository: Distributed via site / CPC Library
- License: GPL (General Public License) or similar open academic license

## Overview
KSSOLV (Kohn-Sham SOLVer) is a MATLAB toolbox designed for solving the Kohn-Sham equations for electronic structure calculations. It utilizes a plane-wave basis set and pseudopotentials. While originally developed for prototyping and educational purposes, KSSOLV 2.0 has evolved into a capable tool for research, enabling algorithm development and direct comparison with standard codes like Quantum ESPRESSO.

**Scientific domain**: Algorithm development, educational DFT, electronic structure
**Target user community**: Developers of DFT algorithms, students learning DFT, researchers needing a flexible prototyping environment

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets
- Norm-conserving pseudopotentials
- Kohn-Sham formulation
- Self-Consistent Field (SCF) methods
- Linear response TDDFT
- GW approximation (in advanced modules)

## Capabilities
- Ground-state energy and density
- Geometry optimization and atomic relaxation
- Ab initio Molecular Dynamics (AIMD)
- Time-Dependent DFT (TDDFT) - real-time and linear response
- Phonon calculations
- GPU acceleration (KSSOLV-GPU)
- Visualization of wavefunctions and densities (via MATLAB)

## Key Strengths

### Ease of Development:
- Written entirely in MATLAB
- Object-oriented design makes code readable and modifiable
- Ideal for testing new SCF mixers, preconditioners, or eigensolvers

### Accessibility:
- leverages MATLAB's robust linear algebra libraries
- Immediate visualization and debugging capabilities within MATLAB environment

### Modern Features:
- Version 2.0 includes significant performance and feature upgrades
- GPU support for accelerated calculations

## Inputs & Outputs
- **Input formats**:
  - MATLAB scripts/functions defining the system
  - Pseudopotential files
  
- **Output data types**:
  - MATLAB data structures (energies, forces, density arrays)
  - Plots and visualizations directly in MATLAB figures

## Interfaces & Ecosystem
- **MATLAB**: Fully integrated into the MATLAB ecosystem.
- **Comparison**: Can use compatible pseudopotentials to compare results with compiled codes.

## Computational Cost
- **Standard**: Slower than compiled codes (C/Fortran) due to MATLAB overhead, though core linear algebra is optimized (BLAS/LAPACK).
- **GPU**: KSSOLV 2.0 GPU implementation offers significant speedups, making it viable for medium-sized research problems.
- **Memory**: Higher memory footprint than optimized Fortran codes.

## Best Practices

### Prototyping Workflow:
- **Small Systems**: Develop and test new algorithms on small molecular systems (e.g., silane, benzene) before scaling up.
- **Profiling**: Use MATLAB's built-in Profiler (`profile on`, `profile viewer`) to identify bottlenecks in new implementations.

### Performance:
- **Vectorization**: When extending KSSOLV, ensure custom code is fully vectorized to maintain MATLAB performance.
- **GPU Usage**: If available, utilize the `KSSOLV-GPU` branch for production-grade SCF cycles.

## Community and Support
- **Hosting**: Distributed via academic sites and CPC Library.
- **Support**: Limited to academic correspondence with authors (C. Yang group).
- **Documentation**: Includes PDFs and inline MATLAB help comments.

## Performance Characteristics
- **Speed**: Generally slower than optimized Fortran/C++ codes (VASP, QE) for production runs, but KSSOLV-GPU offers significant speedups.
- **Efficiency**: Very high "developer efficiency" for implementing new ideas.
- **Scaling**: Limited compared to massive MPI codes, but sufficient for small-to-medium systems and method development.

## Limitations & Known Constraints
- **Language**: Requires MATLAB license (proprietary).
- **Performance**: Not intended for massive high-throughput production runs on supercomputers.
- **Scope**: Feature set is smaller than mature packages like VASP/QE.

## Comparison with Other Codes
- **vs Quantum ESPRESSO**: KSSOLV is for prototyping algorithms that might later be implemented in QE. QE is for production.
- **vs DFT++**: Both emphasize modularity/education, but KSSOLV uses MATLAB while DFT++ uses C++.

## Verification & Sources
**Primary sources**:
1. Official Website: http://kssolv.org/
2. Yang et al., "KSSOLV 2.0..." (arXiv/Journal publications)
3. "KSSOLV - a MATLAB toolbox for solving the Kohn-Sham equations"

**Confidence**: CONFIRMED - Well-documented academic software.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: DFT/MATLAB
- Key Feature: Prototyping/Education
