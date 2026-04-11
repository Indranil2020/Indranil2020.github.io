# aiida-bigdft

## Official Resources
- Source Repository: https://github.com/BigDFT-group/aiida-bigdft-plugin-legacy
- Documentation: Included in repository
- License: Open source

## Overview
**aiida-bigdft** is an AiiDA plugin for the BigDFT wavelet-based DFT code. It enables running BigDFT calculations within the AiiDA framework with provenance tracking, input management via PyBigDFT, and output parsing.

**Scientific domain**: AiiDA plugin for BigDFT wavelet DFT  
**Target user community**: Researchers using BigDFT with AiiDA workflow management

## Theoretical Methods
- BigDFT input generation via PyBigDFT
- Wavelet parameter management
- Output parsing
- AiiDA workflow integration
- Provenance tracking

## Capabilities (CRITICAL)
- BigDFT calculation submission via AiiDA
- Input via PyBigDFT tools
- Output parsing (energies, forces, structure)
- Workflow automation
- Provenance tracking

**Sources**: GitHub repository

## Key Strengths

### BigDFT Integration:
- PyBigDFT input generation
- Wavelet parameter handling
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management

## Inputs & Outputs
- **Input formats**:
  - BigDFT input parameters (via PyBigDFT)
  - Structure data
  
- **Output data types**:
  - Parsed BigDFT output (via LogFile)
  - Energies, forces, structure
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **BigDFT**: Wavelet DFT code
- **PyBigDFT**: Input generation
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: BigDFT-level
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **BigDFT calculations**: Hours (separate)

## Limitations & Known Constraints
- **BigDFT only**: No other code support
- **AiiDA required**: Must have AiiDA
- **Legacy plugin**: May need updates
- **Limited documentation**: Research code

## Comparison with Other Codes
- **vs aiida-vasp**: Different DFT code, same framework
- **vs BigDFT directly**: aiida-bigdft adds provenance
- **Unique strength**: AiiDA plugin for BigDFT wavelet DFT with PyBigDFT integration and provenance tracking

## Application Areas

### BigDFT Workflows:
- Automated BigDFT calculations
- Wavelet-based DFT workflows
- Linear scaling calculations

### Mixed-Code:
- BigDFT + other codes via AiiDA
- Cross-code validation

## Best Practices

### Setup:
- Install AiiDA and configure BigDFT
- Install PyBigDFT
- Test with simple calculation

## Community and Support
- Open source on GitHub
- BigDFT group maintained
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/BigDFT-group/aiida-bigdft-plugin-legacy

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: AiiDA plugin for BigDFT wavelet DFT with PyBigDFT integration and provenance tracking
