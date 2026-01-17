# CoreProjectedHybrids

## Official Resources
- Homepage: https://github.com/bjanesko/CoreProjectedHybrids
- Source Repository: https://github.com/bjanesko/CoreProjectedHybrids
- License: Open Source
- Author: Prof. Benjamin Janesko (Texas Christian University)

## Overview
CoreProjectedHybrids is a PySCF extension module for performing self-consistent and linear-response TDDFT calculations using core-projected hybrid exchange-correlation functionals. These functionals modify the treatment of exact exchange in the core region, offering improved balance between core and valence electron correlation. The package provides specialized functionals for accurate excited state calculations in systems where core-valence interactions are important.

**Scientific domain**: Hybrid functionals, core-level spectroscopy, LR-TDDFT  
**Target user community**: Researchers developing and testing hybrid XC functionals, core-level spectroscopy

## Theoretical Methods
- Linear-Response TDDFT (LR-TDDFT)
- Self-Consistent Field (SCF) with custom functionals
- Core-Projected Hybrid functionals
- Exact exchange modifications
- Core-valence separation
- PySCF integration

## Capabilities
- Self-consistent DFT with core-projected hybrids
- Linear-response TDDFT excited states
- Core-level excitation calculations
- Valence excitation calculations
- Custom XC functional implementation
- PySCF workflow integration

## Key Strengths

### Novel Functional Approach:
- Core-projected exact exchange
- Improved core-valence balance
- Specialized for problematic systems
- Published methodology

### PySCF Integration:
- Seamless PySCF workflow
- Standard input/output formats
- Compatible with PySCF tools
- Full basis set flexibility

### Dual Capability:
- Ground-state SCF
- Linear-response excited states
- Consistent functional treatment
- Unified framework

## Inputs & Outputs
- **Input formats**:
  - PySCF molecule objects
  - Standard basis set specifications
  - Functional parameters
  
- **Output data types**:
  - Ground-state energies
  - Excitation energies
  - Oscillator strengths
  - Orbital properties

## Interfaces & Ecosystem
- **Core dependency**:
  - PySCF (required)

- **Extends**:
  - PySCF SCF module
  - PySCF TDDFT module

## Performance Characteristics
- **Speed**: Standard PySCF performance
- **Accuracy**: Improved for core-valence problems
- **System size**: Standard molecular DFT limits
- **Memory**: Standard PySCF requirements

## Limitations & Known Constraints
- **Specialized use**: Core-projected hybrids are specialized
- **Documentation**: Research-level documentation
- **Validation**: User should validate for their systems
- **Dependencies**: Requires PySCF

## Comparison with Other Codes
- **vs standard hybrids**: Improved core-valence treatment
- **vs range-separated**: Different approach to exchange
- **vs core-valence separated**: Unified treatment
- **Unique strength**: Core-projected hybrid implementation

## Application Areas

### Core-Level Spectroscopy:
- XAS calculations
- Core excitations
- Core-valence correlation

### Functional Development:
- Testing core-projected hybrids
- Benchmarking against standard hybrids
- Method development

## Best Practices
- Compare with standard hybrids
- Validate on known systems
- Use appropriate basis sets for core
- Document functional parameters

## Community and Support
- Open-source on GitHub
- Academic development (TCU)
- Research-focused tool

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/bjanesko/CoreProjectedHybrids
2. Prof. B.G. Janesko publications

**Confidence**: VERIFIED - Active GitHub repository

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Author: Prof. Benjamin Janesko (TCU)
- Purpose: Core-projected hybrid DFT/TDDFT
- Platform: PySCF extension
