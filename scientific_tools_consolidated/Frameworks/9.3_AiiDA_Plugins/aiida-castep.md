# aiida-castep

## Official Resources
- Source Repository: https://github.com/aiidateam/aiida-common-workflows (CASTEP via common workflows)
- Documentation: https://aiida-common-workflows.readthedocs.io/
- License: Open source

## Overview
**aiida-castep** is an AiiDA plugin for the CASTEP DFT code. It enables running CASTEP calculations within the AiiDA framework, supporting the common workflow interface for standardized relaxation, bands, and EOS calculations.

**Scientific domain**: AiiDA plugin for CASTEP plane-wave DFT  
**Target user community**: Researchers using CASTEP with AiiDA workflow management

## Theoretical Methods
- CASTEP input generation
- Pseudopotential handling (NCP, OTF ultrasoft)
- Output parsing
- AiiDA workflow integration
- Common workflow support

## Capabilities (CRITICAL)
- CASTEP calculation submission via AiiDA
- Common relax workflow support
- Pseudopotential management
- Output parsing (energies, forces, structure)
- Provenance tracking

**Sources**: AiiDA common workflows

## Key Strengths

### CASTEP Integration:
- Plane-wave DFT
- Ultrasoft pseudopotentials
- Common workflow interface
- Provenance tracking

### Common Workflows:
- Standardized relaxation
- Cross-code comparison
- Unified I/O

## Inputs & Outputs
- **Input formats**: CASTEP parameters, structure, pseudopotentials
- **Output data types**: Parsed CASTEP output, provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Framework
- **CASTEP**: DFT code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management
- **Accuracy**: CASTEP-level
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **CASTEP calculations**: Hours (separate)

## Limitations & Known Constraints
- **CASTEP only**: No other code support
- **AiiDA required**: Must have AiiDA
- **CASTEP license**: Commercial code
- **Limited documentation**

## Comparison with Other Codes
- **vs aiida-vasp**: Different DFT code, same framework
- **vs CASTEP directly**: aiida-castep adds provenance
- **Unique strength**: AiiDA plugin for CASTEP with common workflow interface

## Application Areas

### CASTEP Workflows:
- Automated CASTEP calculations
- High-throughput with CASTEP
- Common workflow comparison

### Mixed-Code:
- CASTEP + other codes via common workflows
- Cross-code validation

## Best Practices

### Setup:
- Install AiiDA and configure CASTEP
- Set up pseudopotential data
- Test with simple calculation

## Community and Support
- Open source
- AiiDA community
- Limited documentation

## Verification & Sources
**Primary sources**:
1. AiiDA common workflows: https://github.com/aiidateam/aiida-common-workflows

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (via common workflows)
- Specialized strength: AiiDA plugin for CASTEP with common workflow interface
