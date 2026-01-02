# exciting

## Official Resources
- Homepage: https://exciting-code.org/
- Documentation: https://exciting-code.org/ref/documentation
- Source Repository: https://github.com/exciting/exciting
- License: GNU General Public License v3.0

## Overview
exciting is an all-electron full-potential linearized augmented planewave (FP-LAPW) code for DFT and beyond, with particular strength in optical and excited-state properties. It provides advanced capabilities for TDDFT, GW, and BSE calculations with a modern, open-source codebase.

**Scientific domain**: Optical properties, excited states, spectroscopy, all-electron calculations  
**Target user community**: Researchers studying optical properties, excitations, and spectroscopy

## Theoretical Methods
- Density Functional Theory (DFT)
- Full-potential linearized augmented plane wave (FP-LAPW)
- All-electron (no pseudopotentials)
- LDA, GGA, meta-GGA functionals
- Time-Dependent DFT (TDDFT)
- GW approximation (G₀W₀, GW₀, scGW)
- Bethe-Salpeter Equation (BSE)
- DFT+U for correlated systems
- Spin-orbit coupling
- Hybrid functionals
- Random Phase Approximation (RPA)

## Capabilities (CRITICAL)
- Ground-state electronic structure (all-electron)
- Geometry optimization and relaxation
- Total energy, forces, stress tensors
- Band structure and DOS
- Optical properties via TDDFT
- Frequency-dependent dielectric function
- Optical absorption and reflectivity spectra
- GW quasiparticle energies
- BSE for optical excitations including excitonic effects
- X-ray absorption spectroscopy (XAS)
- Electron energy loss spectroscopy (EELS)
- Magneto-optical Kerr effect (MOKE)
- Phonon calculations via linear response
- Elastic constants
- Electric field gradients
- Hyperfine parameters
- Core-level spectroscopy
- Wannier functions
- Berry phase calculations

**Sources**: Official exciting documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - input.xml (XML-based main input)
  - Structure files (various formats)
  - Species files for atomic data
  
- **Output data types**:
  - INFO.OUT (main output)
  - EVALCORE.OUT, EIGVAL.OUT (eigenvalues)
  - Optical spectra files
  - GW output files
  - BSE exciton data
  - Various property-specific outputs

## Interfaces & Ecosystem
- **Framework integrations**:
  - exciting-plus - extended features
  - Wannier90 interface
  - elk2exciting - conversion from Elk
  
- **Visualization**:
  - XCrySDen compatibility
  - exciting-viz tools
  - Standard plotting utilities
  
- **Post-processing**:
  - exciting-optics - optical spectra analysis
  - exciting-xs - excited states analysis
  - Python-based analysis tools

## Limitations & Known Constraints
- **All-electron cost**: Computationally expensive; ~100-200 atom limit for DFT, smaller for GW/BSE
- **GW/BSE expensive**: Very demanding; limited to smaller systems
- **Learning curve**: XML input and LAPW methods require familiarity
- **Memory**: High for all-electron and many-body calculations
- **Parallelization**: MPI and OpenMP but not as scalable as plane-wave codes
- **Documentation**: Good but still developing
- **Community**: Growing but smaller than WIEN2k or Quantum ESPRESSO
- **Installation**: Requires Fortran compiler, libraries (BLAS, LAPACK, FFTW)
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://exciting-code.org/
2. Documentation: https://exciting-code.org/ref/documentation
3. GitHub repository: https://github.com/exciting/exciting
4. A. Gulans et al., J. Phys.: Condens. Matter 26, 363202 (2014) - exciting code

**Secondary sources**:
1. exciting tutorials and workshops
2. Published optical spectra applications
3. GW/BSE benchmark studies
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (mailing list, GitHub issues)
- Academic citations: >500 (main paper)
- Active development: Regular releases, modern codebase
