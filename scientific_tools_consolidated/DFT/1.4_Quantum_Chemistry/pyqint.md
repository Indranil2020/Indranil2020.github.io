# pyqint

## Official Resources
- Homepage: https://github.com/ifilot/pyqint
- Documentation: https://pyqint.readthedocs.io/
- Source Repository: https://github.com/ifilot/pyqint
- License: MIT License

## Overview
pyqint is a Python-based, teaching-oriented implementation of the Hartree-Fock method and molecular integrals. It provides a transparent interface to fundamental electronic structure components, making it excellent for learning and prototyping.

**Scientific domain**: Educational quantum chemistry, molecular integrals  
**Target user community**: Students and educators learning electronic structure theory

## Theoretical Methods
- Restricted Hartree-Fock (RHF)
- Molecular integrals over Gaussian basis functions
- Self-Consistent Field (SCF)
- Geometry optimization
- Mulliken population analysis

## Capabilities (CRITICAL)
- Overlap integrals
- Kinetic energy integrals
- Nuclear attraction integrals
- Two-electron repulsion integrals (ERI)
- Complete SCF procedure
- Geometry optimization
- Population analysis
- Clear Python interface
- Educational focus
- Modular design

## Key Strengths

### Educational Value:
- Clear, readable code
- Step-by-step implementation
- Documentation
- Learning-focused

### Integral Calculations:
- All one-electron integrals
- Two-electron integrals
- Gaussian basis functions
- Contracted GTOs

### SCF Implementation:
- Standard algorithm
- Convergence handling
- Direct and conventional
- Property calculations

### Accessibility:
- Pure Python (with Cython)
- Easy installation
- Minimal dependencies
- Cross-platform

## Inputs & Outputs
- **Input formats**:
  - Python API
  - Molecular coordinates
  - Basis set specifications
  
- **Output data types**:
  - Energies
  - Orbitals
  - Integrals arrays
  - Population data

## Interfaces & Ecosystem
- **NumPy**: Array operations
- **SciPy**: Numerical methods
- **Cython**: Optional acceleration
- **Standard formats**: XYZ files

## Advanced Features

### Integral Engine:
- Obara-Saika scheme
- Recurrence relations
- Contracted Gaussians
- Normalization

### SCF Algorithm:
- Convergence acceleration
- Energy calculation
- Orbital output
- Property evaluation

## Performance Characteristics
- **Speed**: Adequate for teaching
- **Accuracy**: Standard HF accuracy
- **System size**: Small molecules
- **Implementation**: Python/Cython

## Computational Cost
- **Integrals**: O(N^4) ERIs
- **SCF**: Standard HF scaling
- **Typical**: Small molecules for learning

## Limitations & Known Constraints
- **Production**: Not intended for production
- **Methods**: HF only
- **System size**: Small molecules
- **Features**: Basic functionality

## Comparison with Other Codes
- **vs PySCF**: pyqint simpler, educational
- **vs Fermi.jl**: Python vs Julia
- **vs SlowQuant**: Both educational
- **Unique strength**: Clarity over optimization

## Application Areas

### Education:
- Quantum chemistry courses
- Understanding HF
- Integral theory
- Code modification

### Prototyping:
- Testing ideas
- Algorithm development
- Quick implementations

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ifilot/pyqint
2. ReadTheDocs: https://pyqint.readthedocs.io/
3. Szabo & Ostlund implementations

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, MIT)
- Documentation: ReadTheDocs
- Educational focus: Yes
