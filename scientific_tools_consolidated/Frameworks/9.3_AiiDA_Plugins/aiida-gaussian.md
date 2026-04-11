# aiida-gaussian

## Official Resources
- Source Repository: https://github.com/nanotech-empa/aiida-gaussian
- Documentation: Included in repository
- PyPI: https://pypi.org/project/aiida-gaussian/
- License: Open source (MIT)

## Overview
**aiida-gaussian** is an AiiDA plugin for the Gaussian quantum chemistry code. It enables running Gaussian calculations within the AiiDA framework with full provenance tracking, input management, and output parsing for molecular electronic structure calculations.

**Scientific domain**: AiiDA plugin for Gaussian quantum chemistry  
**Target user community**: Researchers using Gaussian with AiiDA workflow management

## Theoretical Methods
- Gaussian input generation and management
- Gaussian output parsing
- AiiDA workflow integration
- Provenance tracking for Gaussian calculations
- Molecular electronic structure

## Capabilities (CRITICAL)
- Gaussian calculation submission via AiiDA
- Input parameter management
- Output parsing (energies, forces, orbitals)
- Workflow automation
- Provenance tracking

**Sources**: GitHub repository

## Key Strengths

### Gaussian Integration:
- All Gaussian calculation types
- Input parameter handling
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management
- Query capabilities

### Molecular Focus:
- Geometry optimization
- Frequency calculations
- TD-DFT
- Multi-step molecular workflows

## Inputs & Outputs
- **Input formats**:
  - Gaussian input parameters
  - Molecular structure
  
- **Output data types**:
  - Parsed Gaussian output
  - Energies, forces, frequencies
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **Gaussian**: Quantum chemistry code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: Gaussian-level
- **System size**: Molecular
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **Gaussian calculations**: Hours (separate)

## Limitations & Known Constraints
- **Gaussian only**: No other code support
- **AiiDA required**: Must have AiiDA
- **Gaussian license**: Commercial code required
- **Molecular systems**: Not for periodic

## Comparison with Other Codes
- **vs aiida-orca**: Different quantum chemistry code
- **vs Gaussian directly**: aiida-gaussian adds provenance
- **Unique strength**: AiiDA plugin for Gaussian with provenance tracking for molecular calculations

## Application Areas

### Molecular Workflows:
- Automated Gaussian calculations
- Multi-step molecular workflows
- High-throughput molecular screening
- Conformer search

### Mixed-Code:
- Gaussian + DFT periodic codes
- QM/MM workflows
- Multi-level calculations

## Best Practices

### Setup:
- Install AiiDA and configure Gaussian
- Set up basis set data
- Test with simple calculation

### Usage:
- Use workchains for complex workflows
- Validate parsed output
- Use AiiDA query for analysis

## Community and Support
- Open source (MIT)
- PyPI installable
- Developed by nanotech-empa
- AiiDA community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/nanotech-empa/aiida-gaussian

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: AiiDA plugin for Gaussian quantum chemistry with provenance tracking
