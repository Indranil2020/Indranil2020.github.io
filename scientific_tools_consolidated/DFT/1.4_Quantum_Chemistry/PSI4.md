# PSI4

## Official Resources
- Homepage: https://psicode.org/
- Documentation: https://psicode.org/psi4manual/master/
- Source Repository: https://github.com/psi4/psi4
- License: GNU Lesser General Public License v3.0

## Overview
PSI4 is an open-source suite of ab initio quantum chemistry programs designed for efficient, high-accuracy simulations of molecular properties. It emphasizes modern software engineering practices, Python integration, and provides state-of-the-art coupled cluster and density functional methods.

**Scientific domain**: Quantum chemistry, molecular properties, method development  
**Target user community**: Quantum chemists, researchers needing accurate molecular calculations

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- LDA, GGA, meta-GGA, hybrid functionals
- Double-hybrid functionals
- MP2, MP3, MP4, CCSD, CCSD(T)
- Orbital-optimized methods (OO-MP2, OO-CCSD)
- Symmetry-adapted perturbation theory (SAPT)
- Algebraic Diagrammatic Construction (ADC)
- Equation-of-motion coupled cluster (EOM-CC)
- Time-Dependent DFT (TDDFT)
- Dispersion corrections (DFT-D3, DFT-D4)
- Solvation models (PCM)
- Density-fitted and Cholesky decomposition methods

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and transition states
- Vibrational frequencies and thermochemistry
- Excited states (TDDFT, EOM-CC, ADC)
- Intermolecular interactions via SAPT
- Molecular properties (dipole, quadrupole, polarizability)
- NMR chemical shifts
- Optical rotation
- Response properties
- Orbital analysis
- Natural bond orbital (NBO) analysis via interface
- Wavefunction analysis
- Composite methods
- Density fitting for efficiency
- GPU acceleration (experimental modules)
- Python API for custom workflows

**Sources**: Official PSI4 documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Python scripts (native interface)
  - Simple input files
  - XYZ coordinate files
  - Z-matrix input
  
- **Output data types**:
  - Detailed output files
  - Energies, gradients, Hessians
  - Molecular orbitals (Molden format)
  - Checkpoint files
  - Wavefunction objects in Python

## Interfaces & Ecosystem
- **Python integration**:
  - Native Python API
  - Psi4NumPy for educational modules
  - Integration with NumPy, SciPy
  - Jupyter notebook support
  
- **External programs**:
  - NBO interface for orbital analysis
  - CFOUR for high-level CC calculations
  - CheMPS2 for DMRG
  - DKH and X2C for relativistic corrections
  
- **Workflow tools**:
  - OptKing for geometry optimization
  - QCEngine for standardized I/O
  - QCArchive infrastructure

## Limitations & Known Constraints
- **Molecular focus**: Not designed for periodic systems
- **System size**: Limited by coupled cluster scaling; ~50-100 atoms for DFT
- **Basis sets**: Gaussian-type only; quality depends on basis
- **Parallelization**: Threading and limited MPI; varies by method
- **Memory**: High-level methods memory-intensive
- **Learning curve**: Moderate; Python knowledge helpful
- **Documentation**: Excellent but assumes quantum chemistry background
- **Platform**: Linux, macOS; Windows support via WSL

## Verification & Sources
**Primary sources**:
1. Official website: https://psicode.org/
2. Documentation: https://psicode.org/psi4manual/master/
3. GitHub repository: https://github.com/psi4/psi4
4. D. G. A. Smith et al., J. Chem. Phys. 152, 184108 (2020) - PSI4 1.4
5. R. M. Parrish et al., J. Chem. Theory Comput. 13, 3185 (2017) - PSI4 1.1

**Secondary sources**:
1. PSI4 tutorials and workshops
2. Psi4NumPy educational modules
3. QCArchive integration documentation
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Very active (forum, GitHub)
- Academic citations: >1,500 (various versions)
- Active development: Regular releases, modern codebase
