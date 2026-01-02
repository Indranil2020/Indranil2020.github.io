# PySCF

## Official Resources
- Homepage: https://pyscf.org/
- Documentation: https://pyscf.org/user.html
- Source Repository: https://github.com/pyscf/pyscf
- License: Apache License 2.0

## Overview
PySCF is a Python-based quantum chemistry package with an emphasis on ab initio methods for molecules and crystals. It provides a simple, lightweight, and efficient platform for developing and testing quantum chemistry methods with excellent scriptability through Python.

**Scientific domain**: Quantum chemistry, method development, molecular and periodic systems  
**Target user community**: Quantum chemists, method developers, researchers needing Python-scriptable quantum chemistry

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- LDA, GGA, meta-GGA, hybrid functionals
- Range-separated hybrids
- MP2, MP3
- Coupled cluster (CCSD, CCSD(T), EOM-CC)
- Multireference methods (CASCI, CASSCF, MRPT)
- FCI (Full Configuration Interaction)
- Algebraic Diagrammatic Construction (ADC)
- GW approximation
- Random Phase Approximation (RPA)
- Periodic boundary conditions (PBC)
- Relativistic methods (X2C, DKH)
- Time-Dependent DFT (TDDFT)

## Capabilities (CRITICAL)
- Ground-state electronic structure (molecules and solids)
- Geometry optimization
- Vibrational frequencies
- Excited states (TDDFT, EOM-CC, ADC)
- GW for quasiparticles
- Periodic systems (Gamma-point and k-points)
- Molecular properties (dipole, quadrupole, polarizability)
- NMR shielding tensors
- Spin-orbit coupling
- Embedding methods (density matrix embedding)
- Multi-reference calculations
- FCI for benchmarking
- Custom method development via Python
- Interface to many external codes
- GPU acceleration (selected modules)

**Sources**: Official PySCF documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Python scripts (native interface)
  - XYZ coordinate files
  - Molecular input via ASE, pymatgen
  - Checkpoint files for restart
  
- **Output data types**:
  - Python objects (energies, wavefunctions)
  - Molecular orbitals
  - Density matrices
  - Property calculations
  - Checkpoint files (.chk)

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - calculator interface
  - pymatgen - structure handling
  - pySCF-forge - extension modules
  - PySCF-NAO - numerical atomic orbitals
  
- **External codes**:
  - DMRG interfaces (Block, CheMPS2)
  - FCIDUMP for external CI codes
  - Molden format for visualization
  
- **Python ecosystem**:
  - NumPy, SciPy for linear algebra
  - Integration with Jupyter notebooks
  - Custom scripting and automation

## Limitations & Known Constraints
- **Python overhead**: Slower than pure Fortran/C++ for routine calculations
- **Memory**: Python memory management can be inefficient for very large calculations
- **Parallelization**: Threading and MPI but not as optimized as commercial codes
- **Basis sets**: Primarily Gaussian-type; limited to available basis set libraries
- **System size**: Molecular focus; periodic features less mature than solid-state codes
- **Documentation**: Good but assumes Python and quantum chemistry knowledge
- **Learning curve**: Requires Python programming skills
- **Platform**: Primarily Linux/macOS; Windows requires care

## Verification & Sources
**Primary sources**:
1. Official website: https://pyscf.org/
2. Documentation: https://pyscf.org/user.html
3. GitHub repository: https://github.com/pyscf/pyscf
4. Q. Sun et al., WIREs Comput. Mol. Sci. 8, e1340 (2018) - PySCF paper
5. Q. Sun et al., J. Chem. Phys. 153, 024109 (2020) - Recent developments

**Secondary sources**:
1. PySCF tutorials and examples
2. Published method development using PySCF
3. Jupyter notebook examples
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Very active (GitHub, Google group)
- Academic citations: >1,000 (main paper)
- Active development: Continuous updates, large contributor base
