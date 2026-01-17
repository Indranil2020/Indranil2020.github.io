# TOMBO

## Official Resources
- Homepage: https://www.tombo.page/
- Documentation: https://www.tombo.page/doc
- Source Repository: Distributed via website request or internal
- License: Contact developers

## Overview
TOMBO (TOhoku Mixed-Basis Orbitals) is an all-electron ab initio simulation code. It utilizes a unique mixed-basis set approach, combining plane waves with atomic orbitals. This allows it to accurately describe both the localized core electron states and the extended valence/conduction states, making it applicable to a wide variety of systems from clusters to periodic crystals.

**Scientific domain**: Clusters, crystals, surfaces, excited states
**Target user community**: Researchers needing high-accuracy all-electron results and excited state properties (GW/BSE)

## Theoretical Methods
- Density Functional Theory (DFT)
- Mixed Basis approach (Plane Waves + Atomic Orbitals)
- All-Electron formulation
- GW Approximation (GWA)
- Bethe-Salpeter Equation (BSE)
- Time-Dependent DFT (TDDFT)
- Geometric Optimization

## Capabilities
- Electronic structure of isolated and periodic systems
- Optical properties via GW-BSE
- Molecular dynamics (MD)
- Thermodynamic properties
- Accurate treatment of core states

## Key Strengths
### Mixed Basis Set
- Efficiently handles both core (localized) and valence (delocalized) electrons
- Bridge between plane-wave codes and LCAO codes
- Reduced computational cost compared to pure plane-wave all-electron methods for some systems

### Excited States
- Robust implementation of GW approximation and Bethe-Salpeter equation
- Accurate prediction of optical absorption spectra

## Inputs & Outputs
- **Input**: Control files specifying basis limits, structure, and calculation type
- **Output**: Energy eigenvalues, wavefunctions, optical spectra, forces

## Computational Cost
- **Scaling**: Standard cubic scaling for DFT, higher for GW/BSE.
- **Efficiency**: Mixed basis reduces the cutoff energy required compared to pure PW codes.

## Interfaces & Ecosystem
- **Availability**:
  - Windows and Mac executables available for free download.
  - Source code available upon request to developers.
- **Input**:
  - Detailed control files for basis and structure definition.

## Advanced Features
- **Time-Dependent GW (TDGW)**:
  - Cutting-edge method for excited state dynamics.
- **LCAO + Plane Wave**:
  - Unique hybrid basis set bridging gap between efficiency and accuracy.

## Community and Support
- **Citation**: Requires citation of Ono et al., Comp. Phys. Comm. 189, 20 (2015).
- **Documentation**: Tutorial textbook and manual available on website.

## Verification & Sources
**Primary sources**:
1. Official Website: https://www.tombo.page/
2. "TOMBO: All-electron mixed-basis ab initio simulation code"

**Confidence**: VERIFIED
**Status**: Active Research Code
