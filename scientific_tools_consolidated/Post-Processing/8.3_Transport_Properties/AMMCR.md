# AMMCR

## Official Resources
- Source Repository: https://github.com/anup12352/AMMCR
- Preprint: https://ar5iv.labs.arxiv.org/html/1907.08005
- License: See repository

## Overview
AMMCR is a code for computing carrier mobility and conductivity using the Rode algorithm, with inputs prepared from first-principles electronic structure calculations. The repository describes an interface targeting VASP output files.

**Scientific domain**: Carrier mobility, electronic conductivity, electron-phonon/impurity scattering models  
**Target user community**: Semiconductor transport modeling from first principles (VASP-based workflows)

## Theoretical Methods
- Semi-classical charge transport modeling
- Rode iterative solution method for the Boltzmann transport equation (as described by the project)

## Capabilities (CRITICAL)
- Mobility calculation
- Electrical conductivity calculation
- VASP-oriented input workflow (as described in repository)

## Inputs & Outputs
- **Input formats**: VASP output files as required by AMMCR (see repository documentation)
- **Output data types**: Mobility and conductivity vs temperature/carrier concentration (workflow-dependent)

## Interfaces & Ecosystem
- Designed to be used with VASP post-processing

## Limitations & Known Constraints
- Requires careful convergence of electronic structure inputs and scattering model parameters.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/anup12352/AMMCR
2. arXiv:1907.08005

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
