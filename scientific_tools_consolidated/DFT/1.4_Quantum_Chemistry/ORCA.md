# ORCA

## Official Resources
- Homepage: https://orcaforum.kofo.mpg.de/
- Documentation: https://www.faccts.de/docs/orca/current/
- Download: https://orcaforum.kofo.mpg.de/app.php/portal
- License: Free for academic use (registration required)

## Overview
ORCA is a modern, general-purpose quantum chemistry program package featuring extensive capabilities in molecular electronic structure calculations. Developed by Frank Neese and coworkers at the Max Planck Institute, ORCA is known for its user-friendly input, comprehensive methods (DFT to high-level coupled cluster), excellent spectroscopic property calculations, and efficient parallelization. It is particularly strong in transition metal chemistry, spectroscopy, and correlated wavefunction methods.

**Scientific domain**: Quantum chemistry, molecular electronic structure, spectroscopy  
**Target user community**: Computational chemists, spectroscopists, transition metal researchers

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (LDA, GGA, meta-GGA, hybrid, double-hybrid)
- Range-separated functionals
- Dispersion corrections (D3, D4)
- Møller-Plesset perturbation theory (MP2, SCS-MP2)
- Coupled Cluster (CCSD, CCSD(T), DLPNO-CCSD(T))
- Multi-reference methods (CASSCF, NEVPT2, MRCI)
- Multireference methods (CASSCF, NEVPT2, MRCI)
- Time-Dependent DFT (TDDFT)
- EOM-CCSD for excited states
- Dispersion corrections (DFT-D3, DFT-D4)
- Solvation models (CPCM, SMD)
- Relativistic methods (ZORA, DKH)

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and transition states
- Vibrational frequencies and thermochemistry
- Reaction pathways (NEB, IRC)
- Spectroscopic properties (UV-Vis, IR, Raman, NMR, EPR, Mössbauer)
- Magnetic properties (g-tensors, hyperfine coupling, ZFS)
- Optical rotation and CD spectra
- X-ray absorption spectroscopy
- Excited state calculations (TDDFT, EOM-CC)
- Multireference calculations for complex systems
- DLPNO methods for large molecules
- QM/MM calculations
- Periodic boundary conditions (limited)
- Automated composite methods
- Thermodynamics and kinetics

**Sources**: Official ORCA documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Simple input file (keyword-based)
  - XYZ coordinate files
  - Flexible input syntax
  
- **Output data types**:
  - Detailed output file
  - Energies, gradients, Hessians
  - Molecular orbitals
  - Spectroscopic data
  - Property calculations

## Interfaces & Ecosystem
- **Visualization**:
  - orca_plot for orbitals and densities
  - Compatible with Avogadro, Chemcraft, GaussView
  
- **Utilities**:
  - orca_mapspc - spectrum processing
  - orca_vib - vibrational analysis
  - orca_comp - composite methods
  
- **QM/MM**:
  - Interface to AMBER, CHARMM force fields
  
- **Workflow integration**:
  - Python wrappers available
  - Compatible with workflow tools

## Limitations & Known Constraints
- **Registration required**: Free for academics but requires registration
- **Not open-source**: Binaries only
- **Periodic systems**: Limited support compared to solid-state codes
- **System size**: Molecular focus; not optimized for extended systems
- **Parallelization**: Excellent but methods vary
- **Platform**: Linux, macOS, Windows binaries
- **Documentation**: Extensive but requires familiarity with quantum chemistry
- **Learning curve**: Moderate; straightforward input but theory knowledge needed

## Verification & Sources
**Primary sources**:
1. Official website: https://orcaforum.kofo.mpg.de/
2. Manual: https://www.faccts.de/docs/orca/6.0/manual/
3. F. Neese, WIREs Comput. Mol. Sci. 12, e1606 (2022) - ORCA 5.0
4. F. Neese, WIREs Comput. Mol. Sci. 8, e1327 (2018) - Software update

**Secondary sources**:
1. ORCA tutorials and workshops
2. Published applications in chemistry
3. Spectroscopy benchmarks
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Software: Free academic binaries available
- Community support: Very active (forum, tutorials)
- Academic citations: >10,000 (various versions)
- Active development: Regular major releases
