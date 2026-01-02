# FHI-aims

## Official Resources
- Homepage: https://fhi-aims.org/
- Documentation: https://fhi-aims-club.gitlab.io/
- Source Repository: Available to users (registration required)
- License: Proprietary (free for academic use)

## Overview
FHI-aims (Fritz Haber Institute ab initio molecular simulations) is an all-electron, full-potential electronic structure code using numeric atom-centered orbitals (NAO). It provides exceptional accuracy and efficiency for molecules and materials, with particular strengths in van der Waals interactions, hybrid functionals, and advanced beyond-DFT methods including GW and RPA.

**Scientific domain**: All-electron DFT, molecules and materials, high-accuracy calculations  
**Target user community**: Researchers requiring accurate all-electron calculations for diverse systems

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Hybrid functionals (PBE0, HSE, B3LYP, etc.)
- Range-separated hybrids
- van der Waals corrections (TS, MBD, dDsC)
- Many-body dispersion (MBD)
- GW approximation (G₀W₀, scGW)
- Time-Dependent DFT (TDDFT)
- GW approximation
- Random Phase Approximation (RPA)
- Second-order Møller-Plesset (MP2)
- Coupled cluster (CCSD, CCSD(T))
- Spin-orbit coupling
- Non-collinear magnetism

## Capabilities (CRITICAL)
- Ground-state electronic structure (all-electron)
- Geometry optimization and transition states
- Total energy, forces, stress tensors
- Molecular dynamics (NVE, NVT, NPT)
- Band structure and DOS
- Optical properties via TDDFT
- GW quasiparticle energies
- RPA correlation energies
- MP2 and coupled cluster calculations
- Accurate van der Waals interactions (MBD method)
- Vibrational spectroscopy (IR, Raman)
- NMR chemical shifts
- Electric polarizabilities and hyperpolarizabilities
- Thermochemistry
- Solvation models (COSMO, MPE)
- Excited states via TDDFT or Delta-SCF
- Dispersion-corrected DFT
- Electron-phonon coupling
- Linear-scaling DFT for large systems

**Sources**: Official FHI-aims documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - control.in (calculation parameters)
  - geometry.in (atomic structure)
  - Numeric basis set definitions
  
- **Output data types**:
  - aims.out (main output with energies, forces)
  - Geometry optimization trajectories
  - DOS and band structure files
  - Molecular orbital outputs
  - Property-specific outputs

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - calculator interface
  - i-PI - path integral MD
  - Phonopy - phonon calculations
  - LibXC - exchange-correlation functionals
  
- **Workflow tools**:
  - AiiDA-FHI-aims (in development)
  - FireWorks integration possible
  
- **Post-processing**:
  - aims-analyzer tools
  - Visualization utilities

## Limitations & Known Constraints
- **Registration required**: Free for academics but requires registration
- **Not fully open-source**: Source available to registered users only
- **Basis sets**: NAO basis sets require convergence testing
- **Memory**: Can be memory-intensive for large basis sets
- **Hybrid functionals**: Computationally expensive; exact exchange costly
- **System size**: Practical limit ~500-1000 atoms for standard DFT; smaller for hybrid/GW
- **k-point sampling**: Primarily for molecules; solids require care
- **Learning curve**: Basis set selection requires experience
- **Parallelization**: Excellent but requires understanding of distribution schemes
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://fhi-aims.org/
2. Documentation: https://fhi-aims.org/documentation
3. V. Blum et al., Comput. Phys. Commun. 180, 2175 (2009) - FHI-aims code paper
4. A. Tkatchenko et al., Phys. Rev. Lett. 108, 236402 (2012) - MBD van der Waals

**Secondary sources**:
1. FHI-aims tutorials and workshops
2. ASE calculator documentation
3. Published benchmark studies
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE (requires registration for full access)
- Source code: Available to registered users
- Community support: Very active (mailing list, workshops)
- Academic citations: >3,000 (main paper)
- Active development: Regular releases, well-maintained
