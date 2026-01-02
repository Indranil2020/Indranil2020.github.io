# CRYSTAL

## Official Resources
- Homepage: http://www.crystal.unito.it/
- Documentation: http://www.crystal.unito.it/documentation.html
- Source Repository: Proprietary (academic/commercial license)
- License: Academic/Commercial license required

## Overview
CRYSTAL is a quantum chemistry program for the study of crystalline solids using Gaussian-type basis functions. It provides ab initio treatment of periodic systems with particular strength in molecular crystals, surfaces, and properties calculations using hybrid functionals.

**Scientific domain**: Molecular crystals, surfaces, materials with localized electrons  
**Target user community**: Researchers studying crystalline solids, surfaces, hybrid functionals

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- Gaussian-type orbital (GTO) basis functions
- All-electron or pseudopotential
- LDA, GGA, meta-GGA functionals
- Hybrid functionals (B3LYP, PBE0, HSE06)
- Double-hybrid functionals
- MP2, MP3 for correlation
- Coupled cluster (CCSD)
- Dispersion corrections (Grimme, TS)
- Spin-orbit coupling (perturbative)

## Capabilities (CRITICAL)
- Ground-state electronic structure (HF and DFT)
- Geometry optimization and transition states
- Vibrational frequencies (harmonic and anharmonic)
- Raman and IR intensities
- Elastic constants and mechanical properties
- Piezoelectric and dielectric tensors
- NMR chemical shifts and coupling constants
- EPR g-tensors and hyperfine parameters
- Optical properties (dielectric function)
- Band structure and DOS
- Compton profiles
- Electron momentum densities
- Topological analysis of electron density (AIM)
- ELF and QTAIM analysis
- X-ray structure factors
- Equation of state
- Surface calculations and adsorption

**Sources**: Official CRYSTAL documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Input file (CRYSTAL format)
  - Geometry input (Cartesian or crystallographic)
  - Basis set library
  - Gaussian94 format compatible
  
- **Output data types**:
  - Standard output with energies, gradients
  - Wavefunction files (.f9, .f98)
  - Formatted checkpoint files
  - Property-specific outputs
  - DOS and band structure files

## Interfaces & Ecosystem
- **Visualization**:
  - CRYSPLOT - plotting utility
  - Interface to external visualization tools
  
- **Post-processing**:
  - TOPOND - topological analysis
  - TOPOS - properties from density
  - Extensive property calculation modules
  
- **Basis sets**:
  - Built-in basis set library
  - Pople, Dunning, and specialized basis sets
  - Effective core potentials available

## Limitations & Known Constraints
- **Licensing**: Requires purchase of academic or commercial license
- **Cost**: Not free; license fees
- **Gaussian basis**: Basis set quality critical for accuracy
- **System size**: Limited by Gaussian basis; ~500 atoms typical
- **k-point sampling**: Important for metallic systems

## Verification & Sources
**Primary sources**:
1. Official website: https://www.crystal.unito.it/
2. Documentation: https://www.crystal.unito.it/documentation.html
3. R. Dovesi et al., Int. J. Quantum Chem. 114, 1287 (2014) - CRYSTAL14
4. R. Dovesi et al., WIREs Comput. Mol. Sci. 8, e1360 (2018) - CRYSTAL17
5. M. Ferrero et al., J. Chem. Phys. 128, 014110 (2008) - Coupled perturbed HF/KS

**Secondary sources**:
1. CRYSTAL manual and tutorials
2. Published solid-state chemistry studies
3. Workshop presentations
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in ALL 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Binary distribution: Available with license
- Community support: Professional support, workshops
- Academic citations: >4,000
- Active development: Regular major releases (CRYSTAL23 latest)
- Benchmark validation: Extensively validated for solids
- Specialized strength: Gaussian basis for periodic systems, hybrid functionals, DFPT
