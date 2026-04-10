# APOST3D

## Official Resources
- GitHub: https://github.com/mgimferrer/APOST3D
- Documentation: https://github.com/mgimferrer/APOST3D/blob/master/DOCUMENTATION.md
- Repository License: GPL-3.0
- Recent publication: APOST-3D: Chemical concepts from wavefunction analysis (J. Chem. Phys. 2024)

## Overview
APOST3D is an open-source Fortran-based code developed at the Universitat de Girona for extracting chemical concepts from wavefunction analysis. It provides a broad set of chemically motivated descriptors, including atomic populations, bond orders, local spin analysis, and effective oxidation states.

**Scientific domain**: Wavefunction analysis, bond orders, local spin, oxidation-state analysis  
**Target user community**: Quantum chemists seeking chemically interpretable descriptors from wavefunction post-processing

## Theoretical Methods
- Atomic population analysis
- Bond order analysis
- Local spin analysis
- Effective oxidation state analysis
- Wavefunction-derived chemical concept extraction

## Capabilities (CRITICAL)
- Open-source wavefunction-analysis code on GitHub
- Calculates atomic populations and bond orders
- Includes local spin analysis and effective oxidation states
- Command-line executable workflow with dedicated documentation
- Uses formatted checkpoint (`.fchk`) and input files as part of the documented workflow

**Sources**: GitHub repository, project documentation, and public publication summary

## Key Strengths

### Broad Chemical-Concept Toolkit:
- Bond orders
- Population analysis
- Local spin descriptors
- Effective oxidation states

### Open and Documented:
- GPL-3.0 licensed
- Public GitHub repository
- Dedicated documentation file and usage instructions

### Modern Public Availability:
- Active public repository
- Clear build instructions
- Recent literature presence

## Inputs & Outputs
- **Input formats**:
  - `.fchk` files
  - APOST3D input files as documented in the repository

- **Output data types**:
  - Atomic populations
  - Bond orders
  - Local spin analysis results
  - Effective oxidation state descriptors

## Workflow and Usage
1. Compile the `apost3d` executable.
2. Prepare `name-input.fchk` and `name-input.inp` in the working directory.
3. Run `apost3d name-input > name-output.apost`.
4. Inspect the resulting chemical descriptors and analysis output.

## Performance Characteristics
- Compiled Fortran workflow
- Suited to detailed descriptor extraction from prepared wavefunction data
- Broader chemical-concept scope than many single-purpose bonding-analysis tools

## Limitations & Known Constraints
- **Input preparation**: Requires compatible formatted checkpoint and input files
- **Command-line workflow**: Less turnkey than some GUI-oriented tools
- **Scope**: Focused on wavefunction-derived descriptors rather than topology visualization

## Comparison with Other Tools
- **vs Chemissian**: APOST3D is a more dedicated scientific descriptor code with open-source workflow and broader chemically motivated analysis scope
- **vs NBO/JANPA**: APOST3D emphasizes a wider set of chemical concepts including local spin and oxidation states in addition to bond orders
- **Unique strength**: Open-source package combining bond orders, local spin, and oxidation-state analysis in one code

## Application Areas
- Bond-order analysis
- Local spin characterization
- Effective oxidation-state studies
- Chemically interpretable wavefunction post-processing

## Community and Support
- Public GitHub repository and issue tracker
- Open-source GPL-3.0 license
- Documented command-line workflow

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mgimferrer/APOST3D
2. Documentation: https://github.com/mgimferrer/APOST3D/blob/master/DOCUMENTATION.md
3. Repository description and usage instructions

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Public repository: ACCESSIBLE
- Documentation: AVAILABLE
- License: GPL-3.0
- Primary use case: Wavefunction analysis for bond orders and related chemical concepts
