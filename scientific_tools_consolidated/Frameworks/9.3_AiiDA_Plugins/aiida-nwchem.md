# aiida-nwchem

## Official Resources
- Source Repository: https://github.com/aiidateam/aiida-nwchem
- Documentation: http://aiida-nwchem.readthedocs.io/
- License: Open source

## Overview
**aiida-nwchem** is an AiiDA plugin for the NWChem quantum chemistry code. It enables running NWChem calculations within the AiiDA framework with provenance tracking, input management, and output parsing for molecular and periodic electronic structure calculations.

**Scientific domain**: AiiDA plugin for NWChem quantum chemistry  
**Target user community**: Researchers using NWChem with AiiDA workflow management

## Theoretical Methods
- NWChem input generation and management
- Basis set handling
- Output parsing
- AiiDA workflow integration
- Provenance tracking

## Capabilities (CRITICAL)
- NWChem calculation submission via AiiDA
- Basis set management
- Output parsing (energies, forces, properties)
- Workflow automation
- Provenance tracking

**Sources**: GitHub repository, ReadTheDocs

## Key Strengths

### NWChem Integration:
- Input parameter handling
- Basis set management
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management

## Inputs & Outputs
- **Input formats**:
  - NWChem input parameters
  - Basis set data
  - Structure data
  
- **Output data types**:
  - Parsed NWChem output
  - Energies, forces, properties
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **NWChem**: Quantum chemistry code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: NWChem-level
- **System size**: Molecular and periodic
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **NWChem calculations**: Hours (separate)

## Limitations & Known Constraints
- **NWChem only**: No other code support
- **AiiDA required**: Must have AiiDA
- **Basis set files**: Need NWChem basis sets
- **Limited documentation**: Research code

## Comparison with Other Codes
- **vs aiida-gaussian**: Different quantum chemistry code
- **vs NWChem directly**: aiida-nwchem adds provenance
- **Unique strength**: AiiDA plugin for NWChem with basis set management and provenance tracking

## Application Areas

### NWChem Workflows:
- Automated NWChem calculations
- High-throughput quantum chemistry
- Molecular property prediction

### Mixed-Code:
- NWChem + DFT codes via AiiDA
- Multi-level calculations

## Best Practices

### Setup:
- Install AiiDA and configure NWChem
- Set up basis set data
- Test with simple calculation

## Community and Support
- Open source on GitHub
- ReadTheDocs documentation
- AiiDA team maintained

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiidateam/aiida-nwchem

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- Specialized strength: AiiDA plugin for NWChem with basis set management and provenance tracking
