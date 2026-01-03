# ComRISB

## Official Resources
- Homepage: https://www.bnl.gov/comscope/software/comscope-software-packages.php
- Documentation: Distributed with ComDMFT
- Source Repository: Part of ComDMFT/Comscope suite
- License: See Comscope project licensing

## Overview
ComRISB (Rotationally Invariant Slave Boson) is a Gutzwiller approximation solver that is part of the Comscope/ComDMFT software suite. It implements the rotationally invariant slave boson method, which provides an alternative to DMFT for treating strong correlations through a variational approach. ComRISB can perform calculations faster than full DMFT while capturing essential correlation effects.

**Scientific domain**: Gutzwiller approximation, rotationally invariant slave boson, strongly correlated materials  
**Target user community**: Researchers studying strongly correlated materials seeking efficient correlation treatments

## Theoretical Methods
- Rotationally Invariant Slave Boson (RISB) method
- Gutzwiller approximation
- Variational approach to correlations
- Alternative to full DMFT calculations
- Mean-field treatment of correlations
- Integration with DFT (DFT+G)

## Capabilities (CRITICAL)
- Gutzwiller-based correlation treatment
- Rotationally invariant formulation
- Faster than full DMFT calculations
- Self-consistent solutions
- Integration with ComDMFT framework
- DFT+Gutzwiller calculations
- Quasiparticle weight calculations
- Orbital occupations and moments
- Works with CyGUTZ implementation

**Sources**: Comscope software packages (https://www.bnl.gov/comscope/), master list notes as part of ComDMFT suite

## Inputs & Outputs
**Input formats**:
- DFT outputs
- Wannier projections
- Interaction parameters
- Configuration files

**Output data types**:
- Quasiparticle weights
- Orbital occupations
- Correlation energies
- Gutzwiller wavefunctions
- Renormalization factors

## Interfaces & Ecosystem
- **ComDMFT**: Distributed as part of ComDMFT
- **CyGUTZ**: Modern Gutzwiller solver implementation
- **Comscope**: Part of BNL Comscope project
- **DFT codes**: Integration via ComDMFT infrastructure

## Limitations & Known Constraints
- Mean-field approximation (not full many-body like DMFT)
- Less accurate than DMFT for some properties
- Does not capture full dynamics
- Zero-temperature formalism primarily
- Limited spectral information compared to DMFT
- Documentation within ComDMFT package
- Requires understanding of Gutzwiller method

## Verification & Sources
**Primary sources**:
1. Comscope website: https://www.bnl.gov/comscope/software/comscope-software-packages.php
2. ComDMFT repository and documentation
3. Master list notes: "RESEARCH CODE - Part of ComDMFT/Comscope suite"

**Secondary sources**:
1. RISB method literature
2. CyGUTZ documentation
3. ComDMFT publications
4. Master list: UNCERTAIN confidence

**Confidence**: UNCERTAIN - Master list marks as research code, part of suite

**Verification status**: ✅ VERIFIED as part of Comscope
- Part of Comscope project: CONFIRMED
- Distributed with ComDMFT: CONFIRMED
- Standalone public repo: NOT FOUND
- Documentation: Within ComDMFT package
- Status: Research code, active within Comscope
