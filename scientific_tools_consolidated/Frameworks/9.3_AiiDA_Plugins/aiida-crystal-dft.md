# aiida-crystal-dft

## Official Resources
- Source Repository: https://github.com/tilde-lab/aiida-crystal-dft
- Documentation: Included in repository
- License: Open source

## Overview
**aiida-crystal-dft** is an AiiDA plugin for the CRYSTAL ab initio code. It enables running CRYSTAL DFT calculations within the AiiDA framework with provenance tracking, input management, and output parsing for periodic and molecular systems.

**Scientific domain**: AiiDA plugin for CRYSTAL DFT calculations  
**Target user community**: Researchers using CRYSTAL with AiiDA workflow management

## Theoretical Methods
- CRYSTAL input generation and management
- Basis set handling
- Output parsing
- AiiDA workflow integration
- Provenance tracking

## Capabilities (CRITICAL)
- CRYSTAL calculation submission via AiiDA
- Basis set management
- Output parsing (energies, forces, frequencies)
- Workflow automation
- Provenance tracking

**Sources**: GitHub repository

## Key Strengths

### CRYSTAL Integration:
- Gaussian-type basis sets
- Input parameter handling
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management

## Inputs & Outputs
- **Input formats**:
  - CRYSTAL input parameters
  - Basis set data
  - Structure data
  
- **Output data types**:
  - Parsed CRYSTAL output
  - Energies, forces, frequencies
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **CRYSTAL**: DFT code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: CRYSTAL-level
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **CRYSTAL calculations**: Hours (separate)

## Limitations & Known Constraints
- **CRYSTAL only**: No other code support
- **AiiDA required**: Must have AiiDA
- **Basis set files**: Need CRYSTAL basis sets
- **Limited documentation**: Research code

## Comparison with Other Codes
- **vs aiida-vasp**: Different DFT code, same framework
- **vs CRYSTAL directly**: aiida-crystal-dft adds provenance
- **Unique strength**: AiiDA plugin for CRYSTAL with Gaussian-type basis set management and provenance tracking

## Application Areas

### CRYSTAL Workflows:
- Automated CRYSTAL calculations
- High-throughput with CRYSTAL
- Frequency calculations

### Mixed-Code:
- CRYSTAL + other codes via AiiDA
- Cross-code validation

## Best Practices

### Setup:
- Install AiiDA and configure CRYSTAL
- Set up basis set data
- Test with simple calculation

## Community and Support
- Open source on GitHub
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/tilde-lab/aiida-crystal-dft

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: AiiDA plugin for CRYSTAL with Gaussian-type basis set management and provenance tracking
