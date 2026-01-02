# Questaal

## Official Resources
- Homepage: https://questaal.org/
- Documentation: https://www.questaal.org/docs/
- Source Repository: https://github.com/questaal/questaal
- License: GNU General Public License v3.0

## Overview
Questaal is a suite of codes for electronic structure calculations using DFT, QSGW (Quasiparticle Self-Consistent GW), and DMFT. It provides both all-electron (LMTO) and pseudopotential implementations with particular strength in strongly correlated systems and GW calculations.

**Scientific domain**: Strongly correlated materials, GW calculations, magnetism  
**Target user community**: Researchers studying correlated electron systems, accurate band structures

## Theoretical Methods
- Density Functional Theory (DFT)
- Linear Muffin-Tin Orbital (LMTO) method
- Full-potential LMTO
- Pseudopotential plane-wave method
- Quasiparticle Self-Consistent GW (QSGW)
- One-shot GW (G₀W₀)
- LDA+DMFT
- LDA, GGA functionals
- DFT+U
- Spin-orbit coupling
- Non-collinear magnetism

## Capabilities (CRITICAL)
- Ground-state electronic structure (all-electron or pseudopotential)
- QSGW for accurate band structures without adjustable parameters
- One-shot GW calculations
- LDA+DMFT for strongly correlated systems
- Geometry optimization
- Total energy and forces
- Band structure and DOS
- Optical properties
- Magnetic properties (moments, exchange interactions)
- Electric field gradients
- Core-level spectroscopy
- Wannier functions
- Transport properties
- Spin dynamics
- Layer Green's function method for surfaces

**Sources**: Official Questaal documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - ctrl file (main control file)
  - site file (structure information)
  - Command-line arguments
  
- **Output data types**:
  - Standard output with energies, moments
  - DOS and band structure files
  - GW self-energy data
  - DMFT output files
  - Property-specific outputs

## Interfaces & Ecosystem
- **Framework integrations**:
  - DMFT solvers interface (CTQMC, etc.)
  - Wannier90 for downfolding
  
- **Visualization**:
  - fplot - plotting utility
  - Standard visualization tools
  
- **Post-processing**:
  - Built-in analysis tools
  - lmf utilities suite

## Limitations & Known Constraints
- **Learning curve**: LMTO methods and input format require familiarity
- **Documentation**: Comprehensive but still evolving
- **Community**: Smaller than major DFT codes
- **QSGW cost**: Computationally expensive; many iterations needed
- **DMFT complexity**: Requires understanding of many-body methods
- **All-electron LMTO**: Limited to ~100-200 atoms typically
- **Parallelization**: MPI supported but scaling varies by method
- **Installation**: Requires libraries (BLAS, LAPACK)
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://questaal.org/
2. Documentation: https://www.questaal.org/docs/
3. GitHub repository: https://github.com/questaal/questaal
4. T. Kotani et al., Phys. Rev. B 76, 165106 (2007) - QSGW method
5. M. van Schilfgaarde et al., Phys. Rev. Lett. 96, 226402 (2006) - QSGW development

**Secondary sources**:
1. Questaal tutorials and examples
2. Published QSGW and DMFT applications
3. Benchmark studies vs experiment
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (tutorials, mailing list)
- Academic citations: >500 (QSGW papers)
- Active development: Regular updates
