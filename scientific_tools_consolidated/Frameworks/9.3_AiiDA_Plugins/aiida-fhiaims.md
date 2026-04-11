# aiida-fhiaims

## Official Resources
- Source Repository: https://github.com/ansobolev/aiida-fhiaims
- Documentation: Included in repository
- License: Open source

## Overview
**aiida-fhiaims** is an AiiDA plugin for the FHI-aims all-electron DFT code. It enables running FHI-aims calculations within the AiiDA framework with provenance tracking, input management, and output parsing.

**Scientific domain**: AiiDA plugin for FHI-aims all-electron DFT  
**Target user community**: Researchers using FHI-aims with AiiDA workflow management

## Theoretical Methods
- FHI-aims input generation and management
- Species defaults handling
- Output parsing
- AiiDA workflow integration
- Provenance tracking

## Capabilities (CRITICAL)
- FHI-aims calculation submission via AiiDA
- Species defaults management
- Output parsing (energies, forces, structure)
- Workflow automation
- Provenance tracking

**Sources**: GitHub repository

## Key Strengths

### FHI-aims Integration:
- Species defaults handling
- Input parameter management
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management

## Inputs & Outputs
- **Input formats**:
  - FHI-aims input parameters
  - Species defaults
  - Structure data
  
- **Output data types**:
  - Parsed FHI-aims output
  - Energies, forces, structure
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **FHI-aims**: All-electron DFT code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: FHI-aims-level
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **FHI-aims calculations**: Hours (separate)

## Limitations & Known Constraints
- **FHI-aims only**: No other code support
- **AiiDA required**: Must have AiiDA
- **Species defaults**: Need FHI-aims species files
- **Limited documentation**: Research code

## Comparison with Other Codes
- **vs aiida-vasp**: Different DFT code, same framework
- **vs FHI-aims directly**: aiida-fhiaims adds provenance
- **Unique strength**: AiiDA plugin for FHI-aims with species defaults management and provenance tracking

## Application Areas

### FHI-aims Workflows:
- Automated FHI-aims calculations
- High-throughput with FHI-aims
- All-electron calculations

### Mixed-Code:
- FHI-aims + other codes via AiiDA
- Cross-code validation

## Best Practices

### Setup:
- Install AiiDA and configure FHI-aims
- Set up species defaults
- Test with simple calculation

## Community and Support
- Open source on GitHub
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ansobolev/aiida-fhiaims

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: AiiDA plugin for FHI-aims with species defaults management and provenance tracking
