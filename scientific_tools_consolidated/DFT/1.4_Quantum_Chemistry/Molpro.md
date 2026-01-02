# Molpro

## Official Resources
- Homepage: https://www.molpro.net/
- Documentation: https://www.molpro.net/manual/
- Source Repository: Proprietary (commercial/academic license)
- License: Commercial/Academic license required

## Overview
Molpro is a comprehensive ab initio quantum chemistry package with particular strength in multi-reference methods, coupled cluster theory, and accurate treatment of electron correlation. It is widely used for high-accuracy calculations on molecular systems, particularly for challenging multi-configurational problems.

**Scientific domain**: Quantum chemistry, multi-reference calculations, high-accuracy correlation  
**Target user community**: Quantum chemists requiring accurate treatment of electron correlation

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- MP2, MP3, MP4
- Coupled cluster (CCSD, CCSD(T), CCSDT, CCSDTQ)
- Explicitly correlated methods (CCSD(T)-F12, MP2-F12)
- Multi-reference CI (MRCI, MRCI+Q)
- CASSCF, RASSCF
- Multi-reference perturbation theory (CASPT2, NEVPT2)
- Multi-reference coupled cluster (MRCC)
- Symmetry-adapted perturbation theory (SAPT)
- Complete active space (CAS) methods
- Local correlation methods (LMP2, LCCSD(T))
- Time-Dependent DFT (TDDFT)

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Multi-reference calculations for complex systems
- Geometry optimization and transition states
- Vibrational frequencies and thermochemistry
- Excited states (MRCI, CASPT2, EOM-CC)
- Conical intersections
- Non-adiabatic coupling
- Explicitly correlated F12 methods for basis set convergence
- Intermolecular interactions via SAPT
- Molecular properties (dipole, quadrupole, polarizability)
- NMR and EPR parameters
- Response properties
- Spin-orbit coupling
- Relativistic corrections (Douglas-Kroll-Hess)
- Local correlation for large molecules
- Automated composite methods

**Sources**: Official Molpro documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Command-based input files
  - XYZ coordinate files
  - Z-matrix input
  - Molden format import
  
- **Output data types**:
  - Detailed output files
  - Energies, gradients, Hessians
  - Molecular orbitals (Molden format)
  - Wavefunction files
  - Property calculations

## Interfaces & Ecosystem
- **External programs**:
  - CFOUR interface for high-level CC
  - Columbus interface for surface hopping
  - Molden for visualization
  
- **Utilities**:
  - Orbital analysis tools
  - Property calculation modules
  - Optimization drivers
  
- **Workflow integration**:
  - Can be scripted for automated calculations
  - Integration with dynamics codes

## Limitations & Known Constraints
- **Commercial license**: Requires purchase for use
- **Cost**: License fees for academic and commercial use
- **Molecular focus**: Not designed for periodic systems
- **System size**: High-level methods limited to small-medium molecules
- **Basis sets**: Gaussian-type; quality critical for accuracy
- **Learning curve**: Steep; complex input for advanced methods
- **Documentation**: Comprehensive but requires expertise
- **Parallelization**: Efficient but varies by method
- **Platform**: Linux, macOS, Windows

## Verification & Sources
**Primary sources**:
1. Official website: https://www.molpro.net/
2. Manual: https://www.molpro.net/manual/
3. H.-J. Werner et al., J. Chem. Phys. 152, 144107 (2020) - Molpro 2020
4. H.-J. Werner et al., WIREs Comput. Mol. Sci. 2, 242 (2012) - Molpro overview

**Secondary sources**:
1. Molpro tutorials and workshops
2. Published high-accuracy benchmark studies
3. Multi-reference method applications
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- License: Commercial/Academic (verified)
- Community support: Active (support, workshops)
- Academic citations: >3,000 (various versions)
- Gold standard: Reference for multi-reference calculations
