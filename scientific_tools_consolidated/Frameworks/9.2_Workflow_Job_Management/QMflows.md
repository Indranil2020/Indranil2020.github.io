# QMflows

## Official Resources
- Source Repository: https://github.com/SCM-NV/qmflows
- Documentation: https://qmflows.readthedocs.io/
- PyPI: https://pypi.org/project/qmflows/
- License: Open source (LGPL-3.0)

## Overview
**QMflows** is a Python library for input generation and task handling in computational chemistry workflows. It provides a high-level interface for running DFT and semi-empirical calculations with ADF, DFTB, ORCA, and CP2K, with automatic job management.

**Scientific domain**: Computational chemistry workflow, multi-code job management  
**Target user community**: Researchers running multi-step quantum chemistry calculations with ADF/DFTB/ORCA/CP2K

## Theoretical Methods
- Multi-code input generation (ADF, DFTB, ORCA, CP2K)
- Workflow orchestration
- Automatic job submission
- Result parsing
- Package management

## Capabilities (CRITICAL)
- ADF, DFTB, ORCA, CP2K input generation
- Workflow definition and execution
- Automatic job management
- Result parsing and storage
- Multi-step calculation chains
- Constrained optimization

**Sources**: GitHub repository, ReadTheDocs

## Key Strengths

### Multi-Code:
- ADF (Amsterdam DFT)
- DFTB (Density Functional Tight Binding)
- ORCA (quantum chemistry)
- CP2K (mixed Gaussian/plane wave)
- Unified API across codes

### Workflow:
- Python-based workflow definition
- Automatic job submission
- Result parsing
- Error handling

### Integration:
- SCM (Software for Chemistry & Materials)
- NAMD (non-adiabatic molecular dynamics)
- Multi-scale workflows

## Inputs & Outputs
- **Input formats**: Molecular structures, calculation parameters
- **Output data types**: Parsed results, energies, gradients, properties

## Interfaces & Ecosystem
- **ADF**: DFT calculations
- **DFTB**: Semi-empirical
- **ORCA**: Quantum chemistry
- **CP2K**: Mixed basis DFT
- **Python**: Core language

## Performance Characteristics
- **Speed**: Workflow management (fast)
- **Accuracy**: Code-dependent
- **System size**: Molecular
- **Automation**: Full

## Computational Cost
- **Framework**: Negligible
- **Calculations**: Hours (separate)

## Limitations & Known Constraints
- **Specific codes**: ADF, DFTB, ORCA, CP2K only
- **Molecular focus**: Primarily molecular systems
- **ADF license**: Commercial code required for ADF
- **Learning curve**: Multi-code setup

## Comparison with Other Codes
- **vs quacc**: QMflows is chemistry-focused, quacc is broader
- **vs AiiDA**: QMflows is lighter, AiiDA has full provenance
- **vs atomate2**: QMflows is molecular, atomate2 is materials
- **Unique strength**: Multi-code computational chemistry workflow with ADF/DFTB/ORCA/CP2K unified API

## Application Areas

### Computational Chemistry:
- Automated DFT calculations
- Multi-step molecular workflows
- Geometry optimization chains
- Spectroscopy calculations

### Multi-Scale:
- DFTB pre-optimization + DFT refinement
- QM/MM workflows
- Conformational search
- Reaction pathway calculation

## Best Practices

### Setup:
- Install supported codes
- Configure job templates
- Test with simple calculations
- Use packages for organization

## Community and Support
- Open source (LGPL-3.0)
- PyPI installable
- SCM maintained
- ReadTheDocs documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/SCM-NV/qmflows

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: Multi-code computational chemistry workflow with ADF/DFTB/ORCA/CP2K unified API
