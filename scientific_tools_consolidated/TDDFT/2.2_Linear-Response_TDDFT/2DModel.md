# 2DModel

## Official Resources
- Homepage: https://github.com/UllrichDFT/2DModel
- Source Repository: https://github.com/UllrichDFT/2DModel
- License: Open Source
- Author: Prof. Carsten Ullrich (University of Missouri)

## Overview
2DModel is a Python code for performing DFT and TDDFT simulations on 2D model solids. Developed by Prof. Carsten Ullrich's group (a leading expert in TDDFT theory), it provides a testbed for exploring fundamental aspects of density functional theory and time-dependent DFT in reduced dimensionality. The code is valuable for methodology development, testing new XC functionals, and understanding TDDFT physics in a controlled setting.

**Scientific domain**: 2D model systems, TDDFT methodology development, linear response theory  
**Target user community**: TDDFT methodology developers, students learning DFT/TDDFT, theoretical physicists

## Theoretical Methods
- Density Functional Theory (DFT)
- Time-Dependent DFT (TDDFT)
- Linear Response TDDFT
- 2D periodic systems
- Model Hamiltonians
- XC functional testing

## Capabilities
- Ground-state DFT for 2D models
- Linear-response TDDFT
- Real-time TDDFT propagation
- Optical response calculations
- Model solid simulations
- Parameter studies
- Methodology benchmarking

## Key Strengths

### Methodology Development:
- Simplified 2D test cases
- Controlled parameter space
- Exact solutions for comparison
- XC functional testing ground

### Expert Authorship:
- Carsten Ullrich (TDDFT textbook author)
- Rigorous theoretical foundation
- Research-grade implementation

### Educational Value:
- Clear Python implementation
- Reduced complexity vs 3D
- Direct theory-to-code correspondence
- Ideal for learning TDDFT

### Flexible Input:
- Configurable parameters
- Multiple calculation modes
- Print output control
- Python-based configuration

## Inputs & Outputs
- **Input formats**:
  - Python input file (infile.py)
  - Configurable parameters:
    - PrintOut: Screen output control
    - Ground state parameters
    - Linear response parameters
  
- **Output data types**:
  - Ground state energies
  - Response functions
  - Optical properties
  - Screen or file output

## Interfaces & Ecosystem
- **Dependencies**:
  - Python standard libraries
  - NumPy/SciPy (typical)

- **Related work**:
  - Ullrich TDDFT textbook
  - Working_TDDFT tutorials

## Theoretical Background
The code implements TDDFT for 2D periodic models, which serve as simplified test systems for understanding:
- Exchange-correlation effects in reduced dimensions
- Linear response formalism
- Collective excitations in 2D
- XC kernel behavior

This follows the pedagogical approach of Prof. Ullrich's textbook "Time-Dependent Density-Functional Theory: Concepts and Applications" (Oxford University Press).

## Performance Characteristics
- **Speed**: Fast (2D models are computationally light)
- **System size**: Model parameters, not atom count
- **Memory**: Minimal requirements
- **Purpose**: Methodology development and teaching

## Limitations & Known Constraints
- **Dimensionality**: 2D only (by design)
- **Real materials**: Model systems, not realistic materials
- **Features**: Focused on fundamental TDDFT
- **Documentation**: Minimal, research-oriented

## Comparison with Other Codes
- **vs production TDDFT codes**: 2DModel for methodology, others for real systems
- **vs 1D model codes**: 2DModel captures additional physics
- **vs full 3D codes**: Much faster, complementary purpose
- **Unique strength**: 2D TDDFT testbed from leading theorist

## Application Areas

### Method Development:
- XC functional benchmarking
- Kernel approximation testing
- Memory effects in TDDFT
- Beyond-adiabatic methods

### Education:
- Teaching TDDFT concepts
- Graduate course projects
- Understanding linear response
- Theory implementation practice

### Fundamental Research:
- 2D material physics insights
- Collective mode studies
- Excitonic effects in 2D
- Correlation in low dimensions

## Best Practices
- Use for understanding, not production calculations
- Compare with analytical limits when available
- Test new XC approximations systematically
- Document parameter choices

## Community and Support
- Open-source on GitHub (UllrichDFT)
- Academic development
- Associated with TDDFT research community
- Python implementation

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/UllrichDFT/2DModel
2. C.A. Ullrich, "Time-Dependent Density-Functional Theory" (OUP, 2012)
3. Associated research publications

**Confidence**: VERIFIED - From established TDDFT research group

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Author: Prof. C.A. Ullrich (leading TDDFT expert)
- Purpose: Methodology development and education
- Language: Python
