# CADPAC

## Official Resources
- Homepage: Historic/Legacy (Cambridge University)
- Documentation: Original manuals and publications
- Note: Cambridge Analytical Derivatives Package
- License: Academic

## Overview
CADPAC (Cambridge Analytic Derivatives Package) is a historic ab initio molecular electronic structure program developed primarily by Nicholas Handy and his group at Cambridge University starting in the 1970s. It pioneered the implementation of analytical first and second derivatives for molecular calculations, fundamentally influencing how modern codes compute gradients and Hessians.

**Scientific domain**: Ab initio quantum chemistry, analytical derivatives  
**Target user community**: Historic significance; development methods now standard in modern codes

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- MP2 perturbation theory with gradients
- DFT (early implementation, one of first)
- Analytical first derivatives (gradients)
- Analytical second derivatives (Hessians)
- Frequency calculations
- CI methods

## Capabilities (CRITICAL)
- First program with analytic second derivatives
- Analytical gradient technology
- Frequency calculations from Hessians
- IR and Raman intensities
- Thermochemical properties
- Force constants
- Normal mode analysis
- Early DFT implementations
- Property calculations
- Influential code design

## Key Strengths

### Historical Significance:
- Pioneered analytical derivatives
- First analytic Hessians in production code
- Influenced Gaussian, GAMESS, others
- Cambridge development (Handy group)
- Seminal papers in derivative theory

### Derivative Technology:
- Full first derivatives
- Full second derivatives
- Efficient algorithms
- Basis set derivatives
- Coupled-perturbed equations

### Properties:
- Vibrational frequencies
- IR intensities
- Raman activities
- Thermochemistry
- Molecular properties

### Software Design:
- Modular structure
- Derivative chain rules
- Efficient coding
- Influenced modern codes

## Inputs & Outputs
- **Input formats**:
  - CADPAC input deck
  - Molecular coordinates
  - Basis set specification
  
- **Output data types**:
  - Energies and gradients
  - Hessian matrices
  - Frequencies
  - Properties

## Interfaces & Ecosystem
- **Standalone**: Self-contained code
- **Output**: Standard formats
- **Integration**: Historic context

## Advanced Features

### Derivative Theory:
- Z-vector technique
- Coupled-perturbed HF
- Gradient optimization
- Hessian-based methods

### Vibrational Analysis:
- Harmonic frequencies
- Normal modes
- Zero-point energies
- Thermodynamic properties

### DFT Development:
- Early functional implementations
- Grid integration
- Exchange-correlation
- Gradient-corrected

## Performance Characteristics
- **Speed**: Competitive for its era
- **Accuracy**: Standard HF/MP2
- **System size**: Small-medium molecules
- **Legacy**: Design principles endure

## Computational Cost
- **HF**: Standard scaling
- **MP2**: O(N^5) gradients
- **Hessians**: Additional overhead
- **Typical**: Research calculations

## Limitations & Known Constraints
- **Legacy code**: No longer developed
- **Availability**: Very limited
- **Platform**: Historic hardware
- **Documentation**: Original papers
- **Superseded**: By modern codes

## Comparison with Other Codes
- **vs Gaussian**: CADPAC influenced Gaussian derivatives
- **vs GAMESS**: Similar influence
- **vs Modern codes**: Derivative algorithms now standard
- **Legacy**: Foundational contribution

## Application Areas

### Historic Applications:
- Geometry optimization
- Transition state location
- Vibrational spectroscopy
- Thermochemistry

### Methodological Impact:
- Analytical derivative algorithm development
- Gradient-corrected DFT
- Coupled-perturbed theory
- Response property theory

## Historical Context

### Development Timeline:
- 1970s: Initial development
- 1980s: Second derivatives
- 1990s: DFT incorporation
- 2000s+: Legacy, methods adopted elsewhere

### Key Publications:
- Pulay force method (independent)
- Handy derivative formulations
- MP2 gradient theory
- DFT implementation papers

## Community and Support
- Historic academic code
- Cambridge University origin
- Handy research group
- Publications as documentation
- Methods now universal

## Verification & Sources
**Primary sources**:
1. Handy et al., original CADPAC papers
2. Cambridge University archives
3. J. Chem. Phys. derivative theory papers
4. Chemical Physics Letters publications

**Confidence**: VERIFIED (Historic)
- Status: Legacy code
- Significance: Pioneer in analytical derivatives
- Impact: Influenced all modern codes
- Methods: Now standard everywhere
