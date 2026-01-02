# Gaussian

## Official Resources
- Homepage: https://gaussian.com/
- Documentation: https://gaussian.com/man/
- Source Repository: Proprietary (commercial license)
- License: Commercial license required

## Overview
Gaussian is one of the most widely-used electronic structure programs in chemistry. It provides a comprehensive suite of methods from Hartree-Fock to high-level coupled cluster with extensive automation and user-friendly interface, making it the de facto standard in computational chemistry.

**Scientific domain**: Computational chemistry, drug design, materials chemistry, spectroscopy  
**Target user community**: Chemists across academia and industry

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- LDA, GGA, meta-GGA, hybrid, double-hybrid functionals
- MP2, MP3, MP4, CCSD, CCSD(T)
- Complete active space (CAS) methods
- Time-Dependent DFT (TDDFT)
- Configuration interaction (CI, QCISD)
- Composite methods (CBS-QB3, G4, W1)
- Solvation models (PCM, SMD, CPCM)
- Dispersion corrections (GD3, GD3BJ)
- ONIOM for QM/MM and multi-layer calculations
- Excited state methods

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and conformational searches
- Transition state searches and IRC calculations
- Vibrational frequencies and thermochemistry
- Excited states (TDDFT, CIS, EOM-CCSD)
- UV-Vis, IR, Raman, NMR, EPR spectra
- Circular dichroism and optical rotation
- Molecular properties (dipole, polarizability, hyperpolarizability)
- Reaction pathways and potential energy surfaces
- ONIOM multi-layer QM/MM calculations
- Solvation effects
- Periodic boundary conditions (limited)
- Automated composite methods for thermochemistry
- Extensive property calculations
- User-friendly input and output

**Sources**: Official Gaussian documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Route section with keywords
  - Z-matrix or Cartesian coordinates
  - Checkpoint files for restart
  - Simple, human-readable format
  
- **Output data types**:
  - Formatted output files
  - Checkpoint files (.chk)
  - Formatted checkpoint files (.fchk)
  - Cube files for densities
  - Archive entries

## Interfaces & Ecosystem
- **Visualization**:
  - GaussView - integrated GUI
  - Molden, Avogadro, Chemcraft compatible
  
- **Workflow integration**:
  - Widely supported by workflow tools
  - Python wrappers available (cclib, GaussianWrangler)
  
- **Analysis tools**:
  - formchk - checkpoint file formatting
  - cubegen - cube file generation
  - freqchk - frequency analysis

## Limitations & Known Constraints
- **Commercial license**: Expensive; requires purchase
- **Cost**: Significant license fees
- **Closed source**: No source code access
- **Molecular focus**: Not optimized for extended systems
- **Periodic systems**: Limited support
- **System size**: Practical limits ~500-1000 atoms for DFT
- **Parallelization**: Efficient but proprietary implementation
- **Platform**: Linux, macOS, Windows (commercial binaries)
- **License management**: Can be restrictive
- **Updates**: Periodic major releases (not continuous)

## Verification & Sources
**Primary sources**:
1. Official website: https://gaussian.com/
2. Manual: https://gaussian.com/man/
3. M. J. Frisch et al., Gaussian 16, Revision C.01, Gaussian, Inc., Wallingford CT, 2016
4. Gaussian White Papers and Technical Notes

**Secondary sources**:
1. Gaussian tutorials and documentation
2. Published applications across all chemistry
3. Textbook references (standard in computational chemistry courses)
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- License: Commercial (verified)
- Community support: Extensive (support, GaussView, tutorials)
- Academic citations: >100,000+ (most cited quantum chemistry code)
- Industry standard: Dominant in pharmaceutical and materials industry
