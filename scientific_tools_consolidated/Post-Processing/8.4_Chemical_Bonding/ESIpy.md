# ESIpy

## Official Resources
- GitHub: https://github.com/jgrebol/ESIpy
- Project Listing: https://quantchemdev.github.io/resources.html
- Documentation/Repository Guidance: GitHub readme and cited Chem. Eur. J. 2024 paper

## Overview
ESIpy is a Python package for calculating population-analysis and aromaticity indicators from different Hilbert-space partitions. It extends the electron-sharing-index ecosystem with a modern Python workflow and can also generate AIMAll-format atomic overlap matrices readable by both ESIpy and the older ESI-3D code.

**Scientific domain**: Electron sharing indices, aromaticity analysis, population analysis  
**Target user community**: Computational chemists studying delocalization, aromaticity, and population analysis in Python workflows

## Theoretical Methods
- Electron sharing indices
- Population analysis from Hilbert-space partitions
- Aromaticity indicators from atomic overlap matrices
- Partition schemes including Mulliken, Löwdin, meta-Löwdin, NAO, and IAO

## Capabilities (CRITICAL)
- Python-based workflow for population and aromaticity analysis
- Supports multiple Hilbert-space partition schemes
- Integrates naturally with PySCF-based calculations
- Can write AIMAll-format AOM directories readable by ESIpy and ESI-3D
- Designed as a modern extension of the ESI/delocalization-analysis ecosystem

**Sources**: GitHub repository, Quantum Chemistry Development Group resources page, and cited project paper

## Key Strengths

### Modern Python Workflow:
- Python-native interface
- PySCF-centered examples
- Easier integration into modern scripting and notebook workflows

### Partition Flexibility:
- Mulliken
- Löwdin and meta-Löwdin
- NAO
- IAO

### Ecosystem Connectivity:
- Bridges newer Python workflows with older ESI-3D-style AOM data
- Useful for delocalization and aromaticity studies
- Supports comparative partition-based analyses

## Inputs & Outputs
- **Input formats**:
  - PySCF molecular and mean-field objects
  - Saved binary molecular-analysis data
  - AIMAll-format AOM directories for compatible workflows

- **Output data types**:
  - Population-analysis results
  - Aromaticity indicators
  - AOM directories for downstream use

## Workflow and Usage
1. Run or load a compatible PySCF calculation.
2. Initialize an `ESI` object with the molecular and calculation data.
3. Choose the partition scheme and target rings or fragments.
4. Print the aromaticity/population results or write AOMs for downstream workflows.

## Performance Characteristics
- Python-oriented workflow suitable for research scripting
- Designed for repeated analysis across different partitions
- Best suited to targeted delocalization and aromaticity studies

## Limitations & Known Constraints
- **Input ecosystem**: Current examples are centered on PySCF-style workflows
- **Method scope**: Focused on population and aromaticity indicators rather than general-purpose wavefunction analysis
- **Specialized audience**: Best for users already interested in ESI-based bonding descriptors

## Comparison with Other Tools
- **vs ESI-3D**: ESIpy provides a more modern Python workflow while preserving interoperability with ESI-3D-style AOMs
- **vs EDDB**: ESIpy emphasizes overlap-matrix-based population and aromaticity indicators, while EDDB emphasizes delocalized density descriptors
- **Unique strength**: Modern Python implementation for electron-sharing and aromaticity analysis across multiple Hilbert-space partitions

## Application Areas
- Aromaticity studies
- Electron delocalization analysis
- Population-analysis comparisons across partition schemes
- Python-based bonding-analysis workflows

## Community and Support
- Public GitHub repository
- Included in the Quantum Chemistry Development Group resources
- Associated with recent literature and active methodological development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/jgrebol/ESIpy
2. Resources page: https://quantchemdev.github.io/resources.html
3. Repository documentation describing PySCF workflows and AOM interoperability

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Public repository: ACCESSIBLE
- Documentation: AVAILABLE
- Primary use case: Python-based electron-sharing and aromaticity analysis
