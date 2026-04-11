# aiida-siesta

## Official Resources
- Source Repository: https://github.com/siesta-project/aiida_siesta_plugin
- Documentation: Included in repository
- License: Open source (MIT)

## Overview
**aiida-siesta** is an AiiDA plugin for the SIESTA DFT code. It enables running SIESTA calculations within the AiiDA framework with provenance tracking, input management, and workflow automation including optical calculations.

**Scientific domain**: AiiDA plugin for SIESTA localized-basis DFT  
**Target user community**: Researchers using SIESTA with AiiDA workflow management

## Theoretical Methods
- SIESTA input generation
- Basis set and pseudopotential handling
- Output parsing
- AiiDA workflow integration
- Optical property calculations

## Capabilities (CRITICAL)
- SIESTA calculation submission via AiiDA
- Basis set and pseudopotential management
- Output parsing (energies, forces, structure)
- Optical calculation support
- Provenance tracking
- WorkChain automation

**Sources**: GitHub repository

## Key Strengths

### SIESTA Integration:
- NAO basis sets
- Pseudopotential handling
- Input parameter management
- Optical calculations
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management

## Inputs & Outputs
- **Input formats**: SIESTA parameters, structure, basis/pseudo
- **Output data types**: Parsed SIESTA output, provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Framework
- **SIESTA**: DFT code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management
- **Accuracy**: SIESTA-level
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **SIESTA calculations**: Hours (separate)

## Limitations & Known Constraints
- **SIESTA only**: No other code support
- **AiiDA required**: Must have AiiDA
- **Basis/pseudo files**: Need SIESTA data
- **Limited documentation**

## Comparison with Other Codes
- **vs aiida-vasp**: Different DFT code, same framework
- **vs SIESTA directly**: aiida-siesta adds provenance
- **Unique strength**: AiiDA plugin for SIESTA with optical calculation support and provenance tracking

## Application Areas

### SIESTA Workflows:
- Automated SIESTA calculations
- High-throughput with SIESTA
- Optical property calculations

### Mixed-Code:
- SIESTA + other codes via AiiDA
- Transiesta transport calculations

## Best Practices

### Setup:
- Install AiiDA and configure SIESTA
- Set up basis and pseudopotential data
- Test with simple calculation

## Community and Support
- Open source (MIT)
- SIESTA project maintained
- INTERSECT project funded

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/siesta-project/aiida_siesta_plugin

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: AiiDA plugin for SIESTA with optical calculation support and provenance tracking
