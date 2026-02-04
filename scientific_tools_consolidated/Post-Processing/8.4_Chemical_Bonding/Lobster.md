# Lobster

## Official Resources
- Homepage: https://www.cohp.de/
- Documentation: http://www.cohp.de/
- Source Repository: Proprietary (free for academic use)
- License: Free for academic research (registration required)

## Overview
Lobster (Local Orbital Basis Suite Towards Electronic-Structure Reconstruction) is a program for chemical bonding analysis based on projecting plane-wave DFT calculations onto a local basis. It enables the calculation of Crystal Orbital Hamilton Populations (COHP), Crystal Orbital Overlap Populations (COOP), and related bonding descriptors from PAW or pseudopotential calculations.

**Scientific domain**: Chemical bonding analysis, electron localization, orbital interactions  
**Target user community**: Chemists and materials scientists studying bonding in solids

## Theoretical Methods
- Projection from plane-wave basis to local atomic orbitals
- Crystal Orbital Hamilton Population (COHP)
- Crystal Orbital Overlap Population (COOP)
- Integrated COHP (ICOHP) for bond strength quantification
- Projected density of states (pDOS)
- Charge analysis (Mulliken, Löwdin)
- Fat band analysis

## Capabilities (CRITICAL)
- COHP analysis for chemical bonding characterization
- COOP analysis for orbital overlap
- ICOHP calculation for quantitative bond strengths
- Projected density of states onto atomic orbitals
- Atom-resolved and orbital-resolved analysis
- Bonding/antibonding character identification
- Charge partitioning (Mulliken, Löwdin populations)
- Fat band structures with orbital character
- Automatic neighbor list generation for bonding pairs
- Interface to VASP (primary)
- Interface to ABINIT
- Interface to Quantum ESPRESSO
- Batch processing capabilities
- Visualization-ready output formats

**Sources**: Official Lobster website, publications, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - POSCAR (VASP structure)
  - POTCAR (PAW information)
  - WAVECAR (VASP wavefunctions)
  - DOSCAR.lobster (optional)
  - lobsterin (Lobster input file)
  
- **Output data types**:
  - COHPCAR.lobster (COHP data)
  - COOPCAR.lobster (COOP data)
  - DOSCAR.lobster (projected DOS)
  - CHARGE.lobster (charge analysis)
  - ICOHPLIST.lobster (integrated COHP values)
  - Orbital-resolved outputs
  - Visualization files

## Interfaces & Ecosystem
- **DFT code interfaces** (verified):
  - VASP - primary and most complete interface
  - Quantum ESPRESSO - supported
  - ABINIT - supported
  
- **Post-processing tools**:
  - LobsterPy - Python package for automated analysis
  - wxDragon - visualization tool for COHP
  - Built-in plotting utilities
  
- **Framework integrations**:
  - pymatgen - can read Lobster outputs
  - Materials Project - uses Lobster for bonding analysis
  - Automated workflows via Python scripts

## Limitations & Known Constraints
- **Projection quality**: Results depend on quality of projection; may not be unique
- **Basis completeness**: Local basis may not span full Hilbert space
- **PAW limitations**: Some approximations in PAW reconstruction
- **Interpretability**: COHP/COOP interpretation requires chemical intuition
- **Code availability**: Free but requires registration; not fully open-source
- **DFT dependency**: Inherits all limitations of underlying DFT calculation
- **Large systems**: Memory and computational cost increase with system size
- **k-point sampling**: Requires adequate k-point mesh from DFT
- **Documentation**: Limited compared to major DFT codes
- **Atomic basis choice**: Results can depend on choice of local basis functions

## Verification & Sources
**Primary sources**:
1. Official website: https://www.cohp.de/
2. S. Maintz et al., J. Comput. Chem. 37, 1030 (2016) - Lobster paper
3. V. L. Deringer et al., J. Phys. Chem. A 115, 5461 (2011) - COHP analysis
4. R. Dronskowski and P. E. Blöchl, J. Phys. Chem. 97, 8617 (1993) - COHP method

**Secondary sources**:
1. VASP wiki on Lobster interface
2. Materials Project documentation on bonding analysis
3. LobsterPy documentation
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE (requires registration for full access)
- Community support: Active (email support, user group)
- Academic citations: >300 (main paper)
- Widely used: Standard tool for bonding analysis in materials science
