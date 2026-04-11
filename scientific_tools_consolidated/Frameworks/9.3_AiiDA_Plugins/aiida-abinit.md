# aiida-abinit

## Official Resources
- Source Repository: https://github.com/aiidateam/aiida-abinit
- Documentation: Included in repository
- License: Open source

## Overview
**aiida-abinit** is an AiiDA plugin for the ABINIT DFT code. It enables running ABINIT calculations within the AiiDA framework with provenance tracking, input management, and output parsing.

**Scientific domain**: AiiDA plugin for ABINIT DFT calculations  
**Target user community**: Researchers using ABINIT with AiiDA workflow management

## Theoretical Methods
- ABINIT input generation and management
- Pseudopotential handling
- Output parsing
- AiiDA workflow integration
- Provenance tracking

## Capabilities (CRITICAL)
- ABINIT calculation submission via AiiDA
- Pseudopotential management
- Output parsing (energies, forces, structure)
- Workflow automation
- Provenance tracking

**Sources**: GitHub repository

## Key Strengths

### ABINIT Integration:
- Input parameter handling
- Pseudopotential management
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management

## Inputs & Outputs
- **Input formats**:
  - ABINIT input parameters
  - Pseudopotential data
  - Structure data
  
- **Output data types**:
  - Parsed ABINIT output
  - Energies, forces, structure
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **ABINIT**: DFT code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: ABINIT-level
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **ABINIT calculations**: Hours (separate)

## Limitations & Known Constraints
- **ABINIT only**: No other code support
- **AiiDA required**: Must have AiiDA
- **Pseudopotential files**: Need ABINIT pseudos
- **Limited documentation**: Research code

## Comparison with Other Codes
- **vs aiida-vasp**: Different DFT code, same framework
- **vs ABINIT directly**: aiida-abinit adds provenance
- **Unique strength**: AiiDA plugin for ABINIT with pseudopotential management and provenance tracking

## Application Areas

### ABINIT Workflows:
- Automated ABINIT calculations
- High-throughput with ABINIT
- DFPT calculations

### Mixed-Code:
- ABINIT + other codes via AiiDA
- Cross-code validation

## Best Practices

### Setup:
- Install AiiDA and configure ABINIT
- Set up pseudopotential data
- Test with simple calculation

## Community and Support
- Open source on GitHub
- AiiDA team maintained
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiidateam/aiida-abinit

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: AiiDA plugin for ABINIT with pseudopotential management and provenance tracking
