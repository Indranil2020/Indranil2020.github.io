# PWDFT (Plane-Wave Density Functional Theory)

## Official Resources
- **Homepage**: https://github.com/ebylaska/PWDFT
- **Source Repository**: https://github.com/ebylaska/PWDFT
- **Developer**: Eric J. Bylaska (Pacific Northwest National Laboratory)
- **License**: Open Source (BSD-style or similar research license implied)

## Overview
PWDFT is a plane-wave density functional theory (DFT) code developed by Eric J. Bylaska at Pacific Northwest National Laboratory (PNNL). It serves as a research platform and mini-application for exploring high-performance computing algorithms, particularly in the context of plane-wave basis sets and pseudopotentials. It is related to the development of NWChem (where Bylaska is a key developer) but exists as a standalone repository for testing and development purposes.

**Scientific domain**: Plane-wave DFT, electronic structure, high-performance computing (HPC)  
**Target user community**: Method developers, HPC researchers, NWChem contributors

## Theoretical Methods
- Kohn-Sham Density Functional Theory
- Plane-wave basis set
- Pseudopotentials (Norm-conserving)
- PSP (Pseudopotential) format support
- Parallel 3D FFTs
- Lagrange multiplier constraints
- Conjugate gradient minimization

## Capabilities
- Ground-state total energy calculations
- Electronic structure minimization
- Wavefunction optimization
- Parallel execution using MPI
- Performance benchmarking for FFTs and parallel data structures
- Mini-app for testing HPC architectures

## Key Strengths
### HPC Research
- Used as a testbed for parallel algorithms
- Investigates scalability of plane-wave methods
- Optimized for memory-bound operations (FFT)

### Connection to NWChem
- Developed by a lead author of NWChem's PSPW module
- Likely shares algorithmic heritage with NWChem's plane-wave implementation

## Inputs & Outputs
- **Input formats**:
  - Input parameter files
  - Pseudopotential files (formatted)
- **Output data types**:
  - Standard output (Energy, convergence)
  - Wavefunction files

## Computational Cost
- **Research Code**: Optimized for benchmarking FFTs and parallel communications.
- **Scaling**: Designed to test limits of HPC scaling; not optimized for production throughput like NWChem.

## Best Practices
- **Usage**: Use for testing new algorithms or compiling on novel architectures (Cray/GPU).
- **Production**: Use NWChem for actual science production runs.

## Comparison with Other Codes
- **vs NWChem**: PWDFT is a smaller, standalone research code/mini-app, while NWChem is a full production suite.
- **vs PWDFT.jl**: PWDFT.jl is a separate educational project in Julia by F. Fathurrahman; this PWDFT is the C++/Fortran research code by Bylaska.

## Verification & Sources
- **Primary Source**: GitHub repository (https://github.com/ebylaska/PWDFT)
- **Developer Info**: Eric Bylaska (PNNL)
- **Confidence**: VERIFIED - Code exists and matches the Master List entry.
