# ALPS

## Official Resources
- Homepage: https://alps.comp-phys.org/
- Documentation: https://alps.comp-phys.org/mediawiki/index.php/Documentation
- Source Repository: https://alps.comp-phys.org/mediawiki/index.php/Download
- License: ALPS License (GPL-compatible, open source)

## Overview
ALPS (Algorithms and Libraries for Physics Simulations) is an international collaboration providing open-source software for simulation of quantum lattice models, including quantum spin systems and strongly correlated electron systems. The project includes both libraries and application codes, with particular relevance for DMFT calculations through its CT-HYB impurity solver implementation.

**Scientific domain**: Quantum lattice models, strongly correlated systems, DMFT  
**Target user community**: Researchers in condensed matter physics studying quantum many-body systems

## Theoretical Methods
- Quantum Monte Carlo algorithms
- Continuous-time quantum Monte Carlo (CT-HYB)
- Exact diagonalization
- Density Matrix Renormalization Group (DMRG)
- Quantum spin models
- Impurity solver for DMFT calculations

## Capabilities (CRITICAL)
- CT-HYB impurity solver for DMFT
- Quantum spin system simulations
- Exact diagonalization for small systems
- DMRG calculations
- C++ libraries with Python interfaces
- Parallel execution support
- HDF5 data format
- Visualization tools
- Integration with DMFT frameworks (DCore, etc.)

**Sources**: Official ALPS website (https://alps.comp-phys.org/), confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- XML parameter files
- Python scripts
- HDF5 archives

**Output data types**:
- Observables and correlation functions
- Green's functions (for CT-HYB solver)
- HDF5 data archives
- Text-based results

## Interfaces & Ecosystem
- **DMFT frameworks**: Used as impurity solver in DCore and other frameworks
- **Python**: Python bindings available
- **Visualization**: Built-in plotting tools
- **Note**: Legacy ALPS project; see ALPSCore for modern continuation

## Limitations & Known Constraints
- Original ALPS project considered legacy
- ALPSCore is the modern continuation with updated libraries
- Installation can be complex
- Documentation scattered across wiki pages
- Some components less actively maintained

## Verification & Sources
**Primary sources**:
1. Official website: https://alps.comp-phys.org/
2. ALPS tutorials and documentation
3. A. Albuquerque et al., J. Magn. Magn. Mater. 310, 1187 (2007) - ALPS project
4. B. Bauer et al., J. Stat. Mech. (2011) P05001 - ALPS updates

**Secondary sources**:
1. ALPS workshops and tutorials
2. Published applications
3. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE (wiki-based)
- Source code: OPEN
- Community: International collaboration
- Status: Legacy project, succeeded by ALPSCore
- CT-HYB solver still used in DMFT applications
