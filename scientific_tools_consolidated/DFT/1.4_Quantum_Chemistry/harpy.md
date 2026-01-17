# harpy

## Official Resources
- Homepage: https://github.com/pwborthwick/harpy
- Documentation: In repository README and code comments
- Source Repository: https://github.com/pwborthwick/harpy
- License: MIT License

## Overview
harpy is a collection of quantum chemistry codes written in Python, focusing on readability and educational value. It includes complete implementations of molecular integrals, Hartree-Fock, and various post-HF methods, designed to help students and educators understand electronic structure theory.

**Scientific domain**: Educational quantum chemistry, molecular integrals, post-HF methods  
**Target user community**: Students, educators, and developers learning quantum chemistry implementations

## Theoretical Methods
- Hartree-Fock (RHF, UHF)
- Molecular integrals (Obara-Saika scheme)
- Møller-Plesset Perturbation Theory (MP2, MP3)
- Coupled Cluster (CCSD, CCSD(T))
- Configuration Interaction (CIS, CISD)
- Time-dependent HF (TDHF/RPA)
- Electron propagator methods
- Various molecular properties

## Capabilities (CRITICAL)
- Complete integral code (Obara-Saika)
- One-electron integrals (overlap, kinetic, nuclear)
- Two-electron repulsion integrals (ERI)
- Direct and conventional SCF
- Post-HF correlation methods
- Response theory implementation
- Readable, well-commented Python code
- Cython acceleration option
- Jupyter notebook examples
- Educational focus throughout

## Key Strengths

### Educational Value:
- Clear, readable implementations
- Extensive code comments
- Step-by-step algorithms
- Learning-focused design
- Easy to modify and extend

### Method Coverage:
- All molecular integrals
- Complete HF implementations
- Popular post-HF methods
- Response properties
- Comprehensive scope

### Code Quality:
- Clean Python style
- Modular design
- Self-documenting
- Cython options for speed
- Jupyter notebook support

### Accessibility:
- Pure Python core
- Minimal dependencies
- Cross-platform
- Easy installation
- Well-organized repository

## Inputs & Outputs
- **Input formats**:
  - Python API
  - XYZ coordinates
  - Basis set specifications
  
- **Output data types**:
  - Energies (HF, MP2, CC)
  - Molecular orbitals
  - Integrals arrays
  - Properties (dipole, etc.)

## Interfaces & Ecosystem
- **NumPy**: Array computations
- **SciPy**: Linear algebra
- **Cython**: Optional acceleration
- **Matplotlib**: Visualization

## Advanced Features

### Integral Engine:
- Obara-Saika recurrence
- Hermite Gaussians
- Boys function
- Schwarz screening
- Contracted GTOs

### Post-HF Methods:
- MP2 and MP3
- CCSD amplitudes
- (T) correction
- EOM-CCSD (in development)
- ADC methods

### Response Theory:
- TDHF/RPA
- CIS excitations
- Transition properties
- Polarizabilities

## Performance Characteristics
- **Speed**: Adequate for teaching
- **Accuracy**: Standard method accuracy
- **System size**: Small molecules
- **Memory**: Standard Python
- **Optimization**: Cython available

## Computational Cost
- **Integrals**: O(N^4) ERIs
- **HF**: Standard SCF scaling
- **Post-HF**: Standard method scaling
- **Typical**: Seconds to minutes for small molecules

## Limitations & Known Constraints
- **Production**: Not for production use
- **System size**: Small molecules only
- **Efficiency**: Educational over optimized
- **Methods**: Standard, no exotic features
- **Community**: Individual project

## Comparison with Other Codes
- **vs PySCF**: harpy simpler, more readable
- **vs pyqint**: Similar educational goals
- **vs SlowQuant**: Different approach
- **vs Production codes**: Educational focus
- **Unique strength**: Clarity and pedagogical value

## Application Areas

### Education:
- Quantum chemistry courses
- Self-study of electronic structure
- Understanding integral theory
- Learning post-HF methods

### Prototyping:
- Testing algorithm ideas
- Quick implementations
- Proof of concept
- Student projects

### Code Development:
- Template for new methods
- Reference implementations
- Debugging aid
- Learning platform

## Best Practices

### Learning:
- Start with integrals
- Progress to SCF
- Then post-HF
- Modify and experiment

### Development:
- Fork and extend
- Add new methods
- Improve efficiency
- Contribute back

## Community and Support
- Open-source MIT license
- GitHub repository
- Individual developer
- Code documentation
- Szabo & Ostlund aligned

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/pwborthwick/harpy
2. Szabo & Ostlund textbook
3. QC integral theory
4. Coupled cluster literature

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, MIT)
- Documentation: README and comments
- Educational value: High
- Implementation: Standard methods
