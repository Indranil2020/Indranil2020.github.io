# NRGLjubljana

## Official Resources
- Homepage: http://nrgljubljana.ijs.si/
- Documentation: http://nrgljubljana.ijs.si/desc.html
- Source Repository: https://github.com/rokzitko/nrgljubljana
- License: GNU General Public License (GPL)

## Overview
NRG Ljubljana is a framework of interrelated computer codes for performing numerical renormalization group (NRG) calculations for quantum impurity problems, such as the Kondo model and Anderson impurity model. It provides highly accurate solutions to impurity problems and can be used as an impurity solver in DMFT calculations, particularly suited for systems with Kondo physics and heavy fermion behavior.

**Scientific domain**: Quantum impurity problems, Kondo physics, DMFT impurity solver  
**Target user community**: Researchers studying quantum impurity problems, heavy fermions, Kondo lattice systems

## Theoretical Methods
- Numerical Renormalization Group (NRG)
- Wilson's NRG approach
- Multi-channel and multi-impurity problems
- Kondo exchange model
- Anderson single and multi-impurity model
- Dynamical quantities (spectral functions)
- Temperature-dependent properties
- Real-frequency calculations

## Capabilities (CRITICAL)
- High-accuracy impurity solver
- Kondo and Anderson impurity models
- Multi-channel impurity problems
- Multi-impurity systems
- Real-frequency spectral functions (no analytical continuation needed)
- Temperature-dependent calculations
- Quantum phase transitions
- Dynamic response functions
- Zero-temperature limit accessible
- DMFT impurity solver capability
- C++ implementation with flexibility

**Sources**: Official NRG Ljubljana website (http://nrgljubljana.ijs.si/), R. Žitko et al., confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- Parameter files for impurity model
- Hybridization functions
- Interaction parameters
- Temperature specifications

**Output data types**:
- Spectral functions (real frequencies)
- Green's functions
- Thermodynamic quantities
- Energy levels and quantum numbers
- Susceptibilities

## Interfaces & Ecosystem
- **DMFT frameworks**: Can be used as impurity solver
- **Standalone**: Independent NRG calculations
- **Research tool**: Precise benchmark calculations

## Limitations & Known Constraints
- Computationally intensive for complex problems
- Requires expertise in NRG methodology
- Setup can be complex
- Memory intensive for large Hilbert spaces
- Best suited for problems with clear energy scales
- Limited to local impurity problems
- Documentation technical

## Verification & Sources
**Primary sources**:
1. Official website: http://nrgljubljana.ijs.si/
2. GitHub repository: https://github.com/rokzitko/nrgljubljana
3. R. Žitko and T. Pruschke, Phys. Rev. B 79, 085106 (2009) - NRG Ljubljana code

**Secondary sources**:
1. NRG Ljubljana documentation and tutorials
2. Published applications in heavy fermion physics
3. Kondo problem benchmarks
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub, GPL)
- Developer: Rok Žitko (Jožef Stefan Institute)
- Well-established NRG implementation
- High-precision impurity solver
