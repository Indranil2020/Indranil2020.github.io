# aiida-common-workflows

## Official Resources
- Source Repository: https://github.com/aiidateam/aiida-common-workflows
- Documentation: https://aiida-common-workflows.readthedocs.io/
- License: Open source (MIT)

## Overview
**aiida-common-workflows** defines a common interface for computational workflows across multiple quantum engines (11+ codes). It provides standardized workflows (relaxation, bands, etc.) that can be run with any supported code through a unified interface.

**Scientific domain**: Common workflow interface across AiiDA plugins, multi-code standardization  
**Target user community**: Researchers wanting code-agnostic AiiDA workflows

## Theoretical Methods
- Common relax workflow (11 quantum engines)
- Common bands workflow
- Common EOS workflow
- Unified input/output specification
- Multi-code interface standardization

## Capabilities (CRITICAL)
- Common relaxation workflow across 11 codes
- Common bands workflow
- Common equation of state workflow
- Standardized input/output
- Multi-code comparison
- Pseudopotential recommendation

**Sources**: GitHub repository

## Key Strengths

### Multi-Code:
- Abinit, BigDFT, CASTEP, CP2K, FLEUR, Gaussian, NWChem, ORCA, QE, Siesta, VASP
- Same workflow, different codes
- Cross-code validation
- Code comparison

### Standardized:
- Common input specification
- Common output format
- Reproducible across codes
- FAIR data principles

### Pseudopotentials:
- SSSP protocol
- PseudoDojo protocol
- Automated recommendation
- Validation

## Inputs & Outputs
- **Input formats**: Standardized structure + parameters
- **Output data types**: Standardized results (energy, forces, structure)

## Interfaces & Ecosystem
- **AiiDA**: Framework
- **11 quantum engines**: Supported codes
- **SSSP/PseudoDojo**: Pseudopotentials
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management
- **Accuracy**: Code-dependent
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Framework**: Negligible
- **Calculations**: Hours (separate)

## Limitations & Known Constraints
- **AiiDA required**: Must have AiiDA
- **Limited workflows**: Relax, bands, EOS
- **Code versions**: Must match supported versions
- **Complex setup**: Multi-code configuration

## Comparison with Other Codes
- **vs individual aiida plugins**: Common interface across all
- **vs quacc**: aiida-common-workflows is AiiDA-native, quacc is multi-engine
- **Unique strength**: Common workflow interface across 11 quantum engines with standardized I/O

## Application Areas

### Cross-Code Validation:
- Compare DFT codes on same system
- Benchmark calculations
- Convergence testing
- Method comparison

### High-Throughput:
- Standardized relaxation
- Automated band structure
- EOS calculations
- Multi-code screening

## Best Practices

### Setup:
- Install AiiDA and required plugins
- Configure pseudopotential families
- Test with simple system
- Compare across codes

## Community and Support
- Open source (MIT)
- AiiDA team maintained
- ReadTheDocs documentation
- Published in Scientific Data

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiidateam/aiida-common-workflows

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Common workflow interface across 11 quantum engines with standardized I/O
