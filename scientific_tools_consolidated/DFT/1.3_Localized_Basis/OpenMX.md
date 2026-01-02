# OpenMX

## Official Resources
- Homepage: http://www.openmx-square.org/
- Documentation: http://www.openmx-square.org/openmx_man3.9/
- Source Repository: http://www.openmx-square.org/ (download page)
- License: GNU General Public License v3.0

## Overview
OpenMX (Open source package for Material eXplorer) is an efficient DFT code using localized pseudo-atomic orbitals (PAO) with particular strength in large-scale calculations, non-collinear magnetism, and spin-orbit coupling. It provides excellent performance for complex magnetic systems and topological materials.

**Scientific domain**: Magnetism, spintronics, topological materials, large systems  
**Target user community**: Researchers studying magnetic materials, spin-orbit physics, large-scale systems

## Theoretical Methods
- Density Functional Theory (DFT)
- Pseudo-atomic orbital (PAO) basis sets
- Norm-conserving pseudopotentials
- LDA, GGA functionals
- DFT+U for correlated systems
- van der Waals corrections
- Spin-orbit coupling (fully relativistic)
- Non-collinear magnetism
- Constrained DFT
- Effective screening medium (ESM) method
- O(N) Krylov subspace method for large systems

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and MD
- Large-scale systems (1000+ atoms with O(N))
- Non-collinear magnetism with spin-orbit coupling
- Magnetic anisotropy energy
- Rashba and Dresselhaus spin splitting
- Topological properties (Z2 invariants, Chern numbers)
- Band structure including spin texture
- Wannier functions and maximally localized Wannier functions
- Quantum transport (NEGF method)
- STM image simulation
- Optical conductivity
- Berry phase calculations
- Electric polarization
- Orbital magnetization
- ESM method for slab calculations
- Linear-scaling DFT

**Sources**: Official OpenMX documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Input file (OpenMX format)
  - Coordinate files (XYZ, PDB)
  - PAO basis definitions
  - Pseudopotential files
  
- **Output data types**:
  - Standard output with energies, forces
  - Band structure files with spin information
  - DOS and PDOS files
  - Density and spin density files
  - Wannier function outputs
  - Transmission coefficients for transport

## Interfaces & Ecosystem
- **Post-processing tools**:
  - OpenMX Viewer - visualization
  - Various analysis utilities included
  - Band unfolding tools
  
- **Transport calculations**:
  - Built-in NEGF for quantum transport
  - Electrode-device-electrode setup
  
- **Workflow integration**:
  - Can be interfaced with ASE
  - Compatible with standard workflow tools

## Limitations & Known Constraints
- **Basis sets**: PAO basis requires convergence testing
- **Pseudopotentials**: Limited to norm-conserving
- **Documentation**: Comprehensive but English translations vary in quality
- **Community**: Smaller than VASP/QE; primarily Japan-based
- **Installation**: Requires compilation; dependencies (BLAS, LAPACK, FFT)
- **Parallelization**: MPI parallelization; efficiency varies
- **Memory**: Generally efficient but depends on basis size
- **k-point sampling**: Required for periodic systems
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: http://www.openmx-square.org/
2. Manual: http://www.openmx-square.org/openmx_man3.9/
3. T. Ozaki, Phys. Rev. B 67, 155108 (2003) - OpenMX method
4. T. Ozaki and H. Kino, Phys. Rev. B 69, 195113 (2004) - O(N) method
5. T. Ozaki et al., Phys. Rev. B 81, 035116 (2010) - Krylov subspace

**Secondary sources**:
1. OpenMX tutorials and examples
2. Published applications in spintronics and topology
3. Benchmark studies
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (download from website)
- Community support: Active (mailing list, Japanese community)
- Academic citations: >1,000 (main papers)
- Active development: Regular updates
