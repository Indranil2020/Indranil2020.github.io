# aiida-openmx

## Official Resources
- Source Repository: https://github.com/azadoks/aiida-openmx
- Documentation: Included in repository
- License: Open source

## Overview
**aiida-openmx** is an AiiDA plugin for the OpenMX DFT code. It enables running OpenMX calculations within the AiiDA framework with provenance tracking, input management using PAO tables, and output parsing.

**Scientific domain**: AiiDA plugin for OpenMX DFT calculations  
**Target user community**: Researchers using OpenMX with AiiDA workflow management

## Theoretical Methods
- OpenMX input generation and management
- PAO (pseudoatomic orbital) table handling
- Output parsing
- AiiDA workflow integration
- Provenance tracking

## Capabilities (CRITICAL)
- OpenMX calculation submission via AiiDA
- PAO table management
- Output parsing (energies, forces, structure)
- Workflow automation
- Provenance tracking

**Sources**: GitHub repository

## Key Strengths

### OpenMX Integration:
- PAO and VPS data management
- Input parameter handling
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management

## Inputs & Outputs
- **Input formats**:
  - OpenMX input parameters
  - PAO/VPS data files
  - Structure data
  
- **Output data types**:
  - Parsed OpenMX output
  - Energies, forces, structure
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **OpenMX**: DFT code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: OpenMX-level
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **OpenMX calculations**: Hours (separate)

## Limitations & Known Constraints
- **OpenMX only**: No other code support
- **AiiDA required**: Must have AiiDA
- **PAO/VPS files**: Need OpenMX data files
- **Limited documentation**: Research code

## Comparison with Other Codes
- **vs aiida-vasp**: Different DFT code, same framework
- **vs OpenMX directly**: aiida-openmx adds provenance
- **Unique strength**: AiiDA plugin for OpenMX with PAO table management and provenance tracking

## Application Areas

### OpenMX Workflows:
- Automated OpenMX calculations
- High-throughput with OpenMX
- Multi-step workflows

### Mixed-Code:
- OpenMX + other codes via AiiDA
- Cross-code validation

## Best Practices

### Setup:
- Install AiiDA and configure OpenMX
- Set up PAO/VPS data
- Test with simple calculation

### Usage:
- Use PAO tables for orbital selection
- Validate parsed output
- Use AiiDA query for analysis

## Community and Support
- Open source on GitHub
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/azadoks/aiida-openmx

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: AiiDA plugin for OpenMX with PAO table management and provenance tracking
