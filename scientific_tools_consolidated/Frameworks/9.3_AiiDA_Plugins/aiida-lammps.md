# aiida-lammps

## Official Resources
- Source Repository: https://github.com/aiidaplugins/aiida-lammps
- Documentation: https://aiida-lammps.readthedocs.io/
- License: Open source (MIT)

## Overview
**aiida-lammps** is an AiiDA plugin for the LAMMPS molecular dynamics code. It enables running LAMMPS calculations within the AiiDA framework with provenance tracking, potential management, and output parsing for classical MD simulations.

**Scientific domain**: AiiDA plugin for LAMMPS molecular dynamics  
**Target user community**: Researchers using LAMMPS with AiiDA workflow management

## Theoretical Methods
- LAMMPS input generation and management
- Potential parameter handling
- Output parsing
- AiiDA workflow integration
- Provenance tracking

## Capabilities (CRITICAL)
- LAMMPS calculation submission via AiiDA
- Potential management (pair_style, pair_coeff)
- Output parsing (energies, forces, trajectory)
- Workflow automation
- Provenance tracking

**Sources**: GitHub repository

## Key Strengths

### LAMMPS Integration:
- Multiple potential types
- Input parameter handling
- Output parsing
- Error recovery

### Provenance:
- Full calculation tracking
- Reproducibility
- Data management

## Inputs & Outputs
- **Input formats**:
  - LAMMPS input parameters
  - Potential data
  - Structure data
  
- **Output data types**:
  - Parsed LAMMPS output
  - Energies, forces, trajectory
  - Provenance graph

## Interfaces & Ecosystem
- **AiiDA**: Workflow framework
- **LAMMPS**: MD code
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: Potential-dependent
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Plugin**: Negligible
- **LAMMPS calculations**: Minutes to hours

## Limitations & Known Constraints
- **LAMMPS only**: No other code support
- **AiiDA required**: Must have AiiDA
- **Potential files**: Need potential data
- **Classical MD**: Not ab initio

## Comparison with Other Codes
- **vs aiida-vasp**: Different code type (MD vs DFT)
- **vs LAMMPS directly**: aiida-lammps adds provenance
- **Unique strength**: AiiDA plugin for LAMMPS with potential management and provenance tracking

## Application Areas

### MD Workflows:
- Automated LAMMPS calculations
- High-throughput MD
- Property prediction from MD

### Mixed-Code:
- DFT + MD via AiiDA
- Multi-scale simulations

## Best Practices

### Setup:
- Install AiiDA and configure LAMMPS
- Set up potential data
- Test with simple calculation

## Community and Support
- Open source (MIT)
- ReadTheDocs documentation
- AiiDA plugins community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiidaplugins/aiida-lammps

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: AiiDA plugin for LAMMPS with potential management and provenance tracking
