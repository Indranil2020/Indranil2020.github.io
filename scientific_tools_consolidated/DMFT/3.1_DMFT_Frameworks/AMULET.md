# AMULET

## Official Resources
- Homepage: http://www.amulet-code.org/
- Documentation: Contact developers for access
- Source Repository: Not publicly available (research code)
- License: Academic/research use

## Overview
AMULET is a DFT+DMFT code developed for studying strongly correlated materials. It represents a research implementation combining density functional theory with dynamical mean field theory. Information about AMULET is limited in public domains, indicating it is primarily a research tool developed and used within specific research groups.

**Scientific domain**: DFT+DMFT, strongly correlated materials  
**Target user community**: Researchers with access to the code through collaborations

## Theoretical Methods
- DFT+DMFT framework
- Dynamical mean field theory
- Electronic structure calculations
- Correlation effects in solids

**Note**: While a research code, documentation is available via the official website.

## Capabilities (CRITICAL)
- **DFT+DMFT Calculations**: Full charge self-consistent implementation
- **Impurity Solvers**:
  - Segment CT-QMC (low temperature efficient)
  - Classical Hirsh-Fye QMC
  - Exact Diagonalization (full rotationally invariant interaction)
- **Disordered Systems**: CPA+DMFT formalism for alloys and chemically disordered compounds
- **Multiple Impurities**: Simultaneous treatment of different correlated shells (e.g., d and f)
- **Magnetic Ordering**: Paramagnetic and any magnetic ordering support
- **Electronic Structure**: Band structure, DOS, magnetic susceptibility
- **Spectroscopy**: ARPES calculation workflow
- **Interfaces**: Quantum Espresso, ELK, TB-LMTO

**Sources**: Official website (http://www.amulet-code.org), MateriApps, zbMATH

## Inputs & Outputs
**Input formats**: Not publicly documented

**Output data types**: Not publicly documented

## Interfaces & Ecosystem
- Integration details not publicly available
- Research code with limited public information


## Performance Characteristics
- **Parallelization**: MPI parallelization
- **Solvers**: Efficient segment solver for low temperatures
- **Efficiency**: CPA implementation for disordered systems

## Comparison with Other Codes
- **vs TRIQS**: AMULET is a specialized DFT+DMFT suite, while TRIQS is a general library
- **vs EDMFTF**: Both handle DMFT, but AMULET has specific strength in CPA+DMFT for alloys
- **vs Wien2k+DMFT**: AMULET integrates with multiple DFT codes (QE, Elk, LMTO)
- **Unique strength**: CPA+DMFT for disordered materials and alloys, ARPES workflow

## Application Areas
- **Strongly Correlated Alloys**: Disordered systems via CPA
- **Magnetic Materials**: Complex magnetic orderings
- **f-electron Systems**: Lanthanides and actinides
- **Surface Science**: ARPES spectra simulations

## Best Practices
- **DFT Interface**: Ensure compatibility with supported DFT codes (QE, Elk)
- **Solver Choice**: Use segment CT-QMC for efficiency where applicable
- **Disorder**: Use CPA module for doped or alloyed systems

## Limitations & Known Constraints
- **Availability**: Primarily a research code, access might be restricted
- **Documentation**: Specialized documentation, may require background knowledge
- **Platform**: Linux environments
- **Updates**: Update frequency lower than major community codes

## Verification & Sources
**Primary sources**:
1. Master list reference (UNCERTAIN confidence level)
2. Limited public information

**Secondary sources**:
1. May exist in research group internal documentation
2. Not confirmed in multiple independent sources

**Confidence**: UNCERTAIN - Master list marks as "UNCERTAIN" confidence

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE (http://www.amulet-code.org)
- Documentation: Available on website
- Source code: Research code (contact developers)
- Capabilities: CPA+DMFT, multiple solvers confirmed
- Citations: Recognized in zbMATH and MateriApps
