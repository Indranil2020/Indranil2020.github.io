# COHP (Crystal Orbital Hamilton Population)

## Official Resources
- Homepage: https://www.cohp.de/
- Documentation: http://www.cohp.de/
- Source Repository: Part of LOBSTER or other implementations
- License: Varies (LOBSTER is academic free)

## Overview
**Status**: METHOD/TOOL - COHP refers to the Crystal Orbital Hamilton Population method, a theoretical framework for analyzing chemical bonding in solids by partitioning the band structure energy into bonding and antibonding contributions. While "COHP" is a method, the term is often used to refer to the LOBSTER code (which implements it) or legacy codes like Cray-COHP.

**Scientific domain**: Chemical bonding analysis, solid-state chemistry  
**Target user community**: Solid-state chemists, materials scientists

## Theoretical Methods
- Partitioning of band structure energy
- Bonding/Antibonding analysis
- Projected density of states (pDOS)
- Crystal Orbital Overlap Population (COOP)
- Integrated COHP (ICOHP) for bond strength

## Capabilities (CRITICAL)
- Visualization of bonding interactions (bonding vs. antibonding)
- Quantitative assessment of bond strengths (ICOHP)
- Analysis of electronic stability
- Implemented in: LOBSTER (standard), tight-binding codes, some DFT packages (SIESTA, LMTO)

**Sources**: Dronskowski and Bloechl, J. Phys. Chem. 97, 8617 (1993)

## Inputs & Outputs
- **Input**: Wavefunctions/Hamiltonians from DFT (VASP, QE, etc.)
- **Output**: COHP vs Energy plots, ICOHP values

## Interfaces & Ecosystem
- **LOBSTER**: The main modern implementation
- **LobsterPy**: Python analysis tool
- **VASP/QE**: Source of electronic structure data

## Application Areas
- Phase stability analysis
- Magnetic ordering origins
- Bond strength in intermetallics
- Surface adsorption

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.cohp.de/
2. Publication: R. Dronskowski, P. E. Bloechl, J. Phys. Chem. 97, 8617 (1993)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Method: STANDARD
- Implementation: Primarily via LOBSTER
- Applications: Bonding analysis
