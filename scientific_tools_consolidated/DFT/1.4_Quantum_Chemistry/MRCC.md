# MRCC

## Official Resources
- Homepage: https://www.mrcc.hu/
- Documentation: https://www.mrcc.hu/manual/
- Source Repository: Available to licensees
- License: Academic and commercial licenses available

## Overview
MRCC (Multi-Reference Correlation Code) is a specialized quantum chemistry program suite featuring arbitrary-order coupled cluster and configuration interaction methods. Developed by Mihály Kállay and collaborators, it is renowned for implementing the highest-order correlation methods available, including automated arbitrary-order CC and CI implementations.

**Scientific domain**: High-order correlation methods, multireference systems, benchmark calculations  
**Target user community**: Researchers requiring very high accuracy or exploring high-order correlation methods

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- Møller-Plesset perturbation theory (MP2-MP∞)
- Coupled Cluster (CC) up to arbitrary order
- CCSD, CCSD(T), CCSDT, CCSDTQ, etc.
- Multireference CC (Mk-MRCC)
- Configuration Interaction (CI) up to full CI
- Linear Response CC
- Equation-of-Motion CC (EOM-CC)
- Symmetry-Adapted Perturbation Theory (SAPT)
- Local correlation methods (LNO-CC)
- Explicitly correlated methods (F12)
- Analytical gradients for many methods
- Relativistic methods (DKH, X2C)
- Automated generation of CC equations

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Very high-order coupled cluster (up to CCSDTQPH)
- Automated coupled cluster code generation
- Geometry optimization with analytic gradients
- Excited states via EOM-CC
- Local correlation for large molecules
- Explicitly correlated F12 methods
- Multi-reference calculations
- Molecular properties
- Benchmark-quality calculations
- Interface to external programs for integral generation
- Efficient implementation of high-level methods
- Parallelization for large calculations

**Sources**: Official MRCC documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - MRCC input file (MINP)
  - XYZ coordinate files
  - Interface inputs from other programs
  
- **Output data types**:
  - Detailed output files
  - Energies, gradients
  - Correlation energies
  - Property calculations
  - Wavefunction analysis

## Interfaces & Ecosystem
- **Integral interfaces**:
  - CFOUR for integral generation
  - ORCA interface
  - Molpro interface
  - Columbus interface
  - Dirac interface for relativistic integrals
  
- **Standalone capabilities**:
  - Can run standalone for many methods
  - Built-in integral code
  
- **Utilities**:
  - dmrcc - main driver
  - Automated equation generation

## Limitations & Known Constraints
- **Registration required**: Free for academics but requires registration
- **Molecular focus**: Not designed for periodic systems
- **System size**: High-order CC limited to very small molecules
- **Memory**: Very high-level methods extremely memory-intensive
- **Computational cost**: High-order CC scales steeply (factorial scaling)
- **Basis sets**: Gaussian-type; large basis required for accuracy
- **Learning curve**: Steep; requires expert knowledge
- **Documentation**: Good but assumes high-level theory knowledge
- **Parallelization**: Efficient but varies by method
- **Platform**: Linux, macOS, Windows

## Verification & Sources
**Primary sources**:
1. Official website: https://www.mrcc.hu/
2. Manual: https://www.mrcc.hu/manual/
3. M. Kállay et al., J. Chem. Phys. 152, 074107 (2020) - MRCC program
4. M. Kállay and P. R. Surján, J. Chem. Phys. 115, 2945 (2001) - High-order CC

**Secondary sources**:
1. MRCC manual and examples
2. Published benchmark calculations
3. High-accuracy thermochemistry studies
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Software: Free for academics (registration required)
- Community support: Active (email support)
- Academic citations: >500
- Unique capability: Automated arbitrary-order CC
