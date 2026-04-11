# aiida-cp2k

## Official Resources
- Source Repository: https://github.com/aiidateam/aiida-cp2k
- Documentation: https://aiida-cp2k.readthedocs.io/
- PyPI: https://pypi.org/project/aiida-cp2k/
- License: Open source (MIT)

## Overview
**aiida-cp2k** is the official AiiDA plugin for the CP2K code. It provides AiiDA-compatible calculation classes, parsers, and workflows for running CP2K calculations within the AiiDA framework, with full provenance tracking and data management.

**Scientific domain**: AiiDA plugin for CP2K DFT/MD calculations  
**Target user community**: Researchers using CP2K with AiiDA workflow management

## Theoretical Methods
- CP2K input generation and management
- CP2K output parsing
- AiiDA workflow integration
- Provenance tracking for CP2K calculations
- Automated error handling

## Capabilities (CRITICAL)
- CP2K calculation submission via AiiDA
- Input parameter management
- Output parsing (energies, forces, structure)
- Workflow automation
- Provenance tracking
- Base workchain for CP2K

**Sources**: GitHub repository, ReadTheDocs

## Key Strengths

### Official Plugin:
- Maintained by AiiDA team
- Regular updates
- Full AiiDA integration
- Standard AiiDA interface

### CP2K Support:
- All CP2K calculation types
- Input parameter handling
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management
- Query capabilities

## Inputs & Outputs
- **Input formats**:
  - CP2K input parameters (dict)
  - Structure data
  - Pseudopotential/basis set data
  
- **Output data types**:
  - Parsed CP2K output
  - Structure data
  - Energy/forces
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **CP2K**: DFT/MD code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: CP2K-level
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **CP2K calculations**: Hours (separate)
- **Typical**: Efficient management

## Limitations & Known Constraints
- **CP2K only**: No other code support
- **AiiDA required**: Must have AiiDA installed
- **CP2K data files**: Need basis/pseudo access
- **Learning curve**: AiiDA + CP2K knowledge needed

## Comparison with Other Codes
- **vs aiida-vasp**: Different DFT code, same framework
- **vs CP2K directly**: aiida-cp2k adds provenance and automation
- **vs aiida-qe**: Different code, similar interface
- **Unique strength**: Official AiiDA plugin for CP2K with full provenance tracking

## Application Areas

### CP2K Workflows:
- Automated CP2K calculations
- High-throughput CP2K
- MD simulations with provenance
- Geometry optimization workflows

### Mixed-Code Workflows:
- CP2K + other codes via AiiDA
- Multi-level calculations
- QM/MM workflows
- Cross-code validation

## Best Practices

### AiiDA Setup:
- Install AiiDA first
- Configure CP2K code in AiiDA
- Set up basis/pseudo data
- Test with simple calculation

### Calculations:
- Use BaseWorkChain for error recovery
- Validate input parameters
- Check parsed output
- Use AiiDA query for analysis

## Community and Support
- Open source (MIT)
- PyPI installable
- ReadTheDocs documentation
- AiiDA team maintained
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiidateam/aiida-cp2k

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- PyPI: AVAILABLE
- Specialized strength: Official AiiDA plugin for CP2K with full provenance tracking
