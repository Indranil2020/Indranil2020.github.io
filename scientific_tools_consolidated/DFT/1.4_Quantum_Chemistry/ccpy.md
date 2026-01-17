# ccpy

## Official Resources
- Homepage: https://github.com/piecuch-group/ccpy
- Documentation: In repository
- Source Repository: https://github.com/piecuch-group/ccpy
- License: GNU General Public License v3.0

## Overview
ccpy is a Python-based coupled-cluster package developed by the Piecuch group at Michigan State University. It implements a variety of ground and excited-state coupled-cluster methods using a hybrid Python-Fortran approach for computational efficiency.

**Scientific domain**: Coupled-cluster theory, excited states, post-HF methods  
**Target user community**: Researchers applying advanced coupled-cluster methods

## Theoretical Methods
- CCD, CCSD, CCSDT, CCSDTQ
- CCSD(T) and other perturbative triples
- Completely Renormalized CC (CR-CC)
- Equation-of-Motion CC (EOM-CC)
- Left eigenstate CC
- Active-space CC methods
- Similarity-transformed EOM

## Capabilities (CRITICAL)
- Ground-state coupled cluster
- Excited-state EOM-CC
- Left eigenstates for properties
- Perturbative corrections
- CR-CC methods (CR-CC(2,3))
- Active-space extensions
- Hybrid Python/Fortran
- Efficient contractions

## Key Strengths

### Method Coverage:
- Multiple CC approximations
- Perturbative corrections
- Renormalized methods
- Excited states

### Piecuch Methods:
- CR-CC(2,3)
- DEA/DIP-EOMCC
- Active-space variants
- Size-extensive corrections

### Implementation:
- Python front-end
- Fortran performance
- Modular design
- Extensibility

### Research Focus:
- Method development
- New approximations
- Benchmark calculations

## Inputs & Outputs
- **Input formats**:
  - Molecular integrals
  - Python scripts
  
- **Output data types**:
  - Energies
  - Amplitudes
  - Excitation energies
  - Transition properties

## Interfaces & Ecosystem
- **PySCF**: Integral interface
- **NumPy**: Array operations
- **Fortran**: Computation kernels

## Advanced Features

### CR-CC Methods:
- Non-iterative corrections
- Size-extensivity
- Multi-reference character
- Bond breaking

### EOM-CC Variants:
- EOM-CCSD
- EOM-CCSDT
- IP/EA/EE variants
- Transition moments

### Active Space:
- CCSDt
- CC(t;3)
- Reduced scaling
- Large systems

## Performance Characteristics
- **Speed**: Fortran-accelerated
- **Accuracy**: High-level methods
- **System size**: Medium molecules
- **Parallelization**: OpenMP

## Computational Cost
- **CCSD(T)**: O(N^7) perturbative
- **EOM-CCSD**: O(N^6) per state
- **CR-CC**: Additional corrections
- **Typical**: Moderate molecules

## Limitations & Known Constraints
- **Documentation**: Research-focused
- **Large systems**: Standard CC limitations
- **User interface**: Requires expertise
- **Community**: Research group centered

## Comparison with Other Codes
- **vs CFOUR**: Both high-level CC; different methods
- **vs MRCC**: Both advanced CC
- **vs ccq**: ccpy more methods, Piecuch focus
- **Unique strength**: CR-CC, Piecuch group methods

## Application Areas

### Excited States:
- EOM-CC calculations
- Transition properties
- Multi-state problems

### Strong Correlation:
- CR-CC for bond breaking
- Active-space methods
- Challenging systems

### Benchmarks:
- Reference calculations
- Method validation
- New developments

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/piecuch-group/ccpy
2. Piecuch group publications
3. Michigan State University

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Academic group: Piecuch group
- Active development: Yes
