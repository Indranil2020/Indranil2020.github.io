# aiida-orca

## Official Resources
- Source Repository: https://github.com/aiidateam/aiida-common-workflows (ORCA via common workflows)
- Documentation: https://aiida-common-workflows.readthedocs.io/
- License: Open source

## Overview
**aiida-orca** is an AiiDA plugin for the ORCA quantum chemistry code. It enables running ORCA calculations within the AiiDA framework, supporting the common workflow interface for standardized calculations.

**Scientific domain**: AiiDA plugin for ORCA quantum chemistry  
**Target user community**: Researchers using ORCA with AiiDA workflow management

## Theoretical Methods
- ORCA input generation
- Basis set handling
- Output parsing
- AiiDA workflow integration
- Common workflow support

## Capabilities (CRITICAL)
- ORCA calculation submission via AiiDA
- Common relax workflow support
- Basis set management
- Output parsing (energies, forces, properties)
- Provenance tracking

**Sources**: AiiDA common workflows

## Key Strengths

### ORCA Integration:
- Multi-reference methods
- TDDFT
- EPR/NMR predictions
- Common workflow interface
- Provenance tracking

### Common Workflows:
- Standardized relaxation
- Cross-code comparison
- Unified I/O

## Inputs & Outputs
- **Input formats**: ORCA parameters, structure, basis sets
- **Output data types**: Parsed ORCA output, provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Framework
- **ORCA**: Quantum chemistry code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management
- **Accuracy**: ORCA-level
- **System size**: Molecular
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **ORCA calculations**: Hours (separate)

## Limitations & Known Constraints
- **ORCA only**: No other code support
- **AiiDA required**: Must have AiiDA
- **ORCA license**: Free for academic use
- **Molecular focus**: Not for periodic

## Comparison with Other Codes
- **vs aiida-gaussian**: Different quantum chemistry code
- **vs ORCA directly**: aiida-orca adds provenance
- **Unique strength**: AiiDA plugin for ORCA with common workflow interface

## Application Areas

### ORCA Workflows:
- Automated ORCA calculations
- High-throughput quantum chemistry
- Spectroscopy calculations

### Mixed-Code:
- ORCA + DFT codes via common workflows
- Multi-level calculations

## Best Practices

### Setup:
- Install AiiDA and configure ORCA
- Set up basis set data
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
- Specialized strength: AiiDA plugin for ORCA quantum chemistry with common workflow interface
