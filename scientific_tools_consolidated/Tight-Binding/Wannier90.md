# Wannier90

## Official Resources
- Homepage: https://wannier.org/
- Documentation: https://wannier.org/user-guide/
- Source Repository: https://github.com/wannier-developers/wannier90
- License: GNU General Public License v2.0

## Overview
Wannier90 is the standard code for computing maximally-localized Wannier functions (MLWFs) from electronic structure calculations. It constructs localized orbitals from Bloch states calculated by DFT codes, enabling efficient calculations of band structures, transport properties, and tight-binding models. It is widely used as a post-processing tool and is integrated with virtually all major DFT packages.

**Scientific domain**: Electronic structure, tight-binding models, localized orbitals, topological materials  
**Target user community**: Researchers needing Wannier functions, tight-binding models, or k-space interpolation

## Theoretical Methods
- Maximally Localized Wannier Functions (MLWFs)
- Marzari-Vanderbilt localization
- Projected Wannier functions
- Disentanglement procedures
- Band structure interpolation
- Berry phase calculations
- Boltzmann transport via Wannier functions
- Orbital magnetization
- Symmetry-adapted Wannier functions

## Capabilities (CRITICAL)
- Generation of maximally localized Wannier functions from DFT
- Wannier interpolation of electronic band structures
- Tight-binding Hamiltonian construction in Wannier basis
- Berry phase calculations (polarization, Chern numbers)
- Berry curvature and anomalous Hall conductivity
- Orbital magnetization
- Gyrotropic effects
- Shift current and nonlinear optical properties
- Interface to DFT+DMFT codes (TRIQS, w2dynamics)
- Interface to electron-phonon codes (EPW, Perturbo)
- Interface to transport codes (BoltzWann)
- Interface to topological analysis tools (WannierTools, Z2Pack)
- Plotting and visualization utilities
- Parallel execution support

**Sources**: Official Wannier90 website, user guide, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - .win file (Wannier90 input file)
  - .mmn file (overlap matrices from DFT)
  - .amn file (projection matrices from DFT)
  - .eig file (eigenvalues from DFT)
  - .dmn file (for Berry curvature calculations)
  - DFT interface outputs (VASP, QE, ABINIT, etc.)
  
- **Output data types**:
  - .wout file (main output with convergence and diagnostics)
  - _hr.dat (Hamiltonian in Wannier basis)
  - _centres.xyz (Wannier function centers)
  - _band.dat, _band.gnu (interpolated band structures)
  - .chk file (checkpoint for restarts)
  - .nnkp file (k-point information for DFT interface)
  - _berry.dat (Berry phase/curvature data)

## Interfaces & Ecosystem
- **DFT code interfaces** (verified):
  - VASP - via VASP2WANNIER90 interface
  - Quantum ESPRESSO - via pw2wannier90.x
  - ABINIT - native support
  - CASTEP - supported
  - FLEUR - Wannier90 interface
  - WIEN2k - via Wien2Wannier interface
  - Siesta - supported
  - FPLO - supported
  - Elk - supported
  - CRYSTAL - supported
  - OpenMX - supported
  - Quantum-ATK - supported
  
- **Downstream applications**:
  - EPW - electron-phonon coupling (Quantum ESPRESSO)
  - PERTURBO - carrier dynamics and transport
  - WannierTools - topological property calculations
  - WannierBerri - Berry phase properties and transport
  - BoltzWann - Boltzmann transport
  - Z2Pack - topological invariants
  - TRIQS/DFTTools - DFT+DMFT downfolding
  - w2dynamics - DMFT calculations
  
- **Framework integrations**:
  - AiiDA - aiida-wannier90 plugin for workflows
  - ASE - structure preparation
  - pymatgen - structure I/O
  
- **Python interfaces**:
  - WannierBerri - Python interface for transport
  - pythtb - TB models from Wannier output
  - TBmodels - Python TB manipulation

## Limitations & Known Constraints
- **Initial guess dependency**: Convergence sensitive to initial projections; requires physical intuition
- **Disentanglement**: Entangled bands require careful disentanglement procedure; can fail for complex band structures
- **Gauge freedom**: Wannier functions not unique; different choices valid but affect localization
- **System size**: Practical k-point mesh limited by DFT calculation size; dense meshes expensive
- **Accuracy**: Interpolation accuracy depends on localization quality and k-point sampling
- **Topological properties**: Requires careful validation; numerical artifacts possible
- **Documentation**: User guide comprehensive but learning curve steep for beginners
- **Memory**: Large systems with many bands can be memory intensive
- **Convergence**: Some systems difficult to converge (metals, partially filled bands)

## Verification & Sources
**Primary sources**:
1. Official website: https://wannier.org/
2. User guide: http://www.wannier.org/support/
3. N. Marzari and D. Vanderbilt, Phys. Rev. B 56, 12847 (1997) - Original MLWF paper
4. I. Souza et al., Phys. Rev. B 65, 035109 (2001) - Wannier interpolation
5. A. A. Mostofi et al., Comput. Phys. Commun. 178, 685 (2008) - Wannier90 code
6. G. Pizzi et al., J. Phys.: Condens. Matter 32, 165902 (2020) - Wannier90 v3

**Secondary sources**:
1. Wannier90 user guide and tutorials
2. Published studies using Wannier90 (>2,500 citations)
3. DFT code documentation
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in ALL 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, GPL v2)
- Community support: Mailing list, GitHub issues
- Academic citations: >3,000
- Active development: Regular releases, active GitHub
- DFT integration: Standard post-processing tool
- Industry standard: Most widely used Wannier function code
- Educational adoption: Standard tool for electronic structure analysis
- Source code: OPEN (GitHub)
- Community support: Active (mailing list, workshops)
- Academic citations: >2,500 (Google Scholar for main papers)
- DFT interfaces: 12+ verified
- Ecosystem tools: Multiple actively maintained packages
