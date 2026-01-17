# GBasis

## Official Resources
- Homepage: https://gbasis.qcdevs.org/
- Documentation: https://gbasis.readthedocs.io/
- Source Repository: https://github.com/theochem/gbasis
- PyPI: https://pypi.org/project/gbasis/
- License: GNU Lesser General Public License v3.0

## Overview
GBasis is a Python library for evaluating and integrating Gaussian-type orbitals and integrals. Part of the QCDevs project and originating from the HORTON project, it provides a modular framework for computing molecular integrals, densities, and related quantities needed for quantum chemistry calculations.

**Scientific domain**: Quantum chemistry, molecular integrals, basis set evaluation  
**Target user community**: Quantum chemistry developers, researchers needing integral evaluation, educational users

## Theoretical Methods
- Gaussian-type orbital (GTO) evaluation
- One-electron integrals
- Two-electron integrals
- Molecular orbital basis transformations
- Density and potential evaluation
- Electrostatic moments

## Capabilities (CRITICAL)
- GTO basis function evaluation
- Overlap integrals
- Kinetic energy integrals
- Nuclear attraction integrals
- Electron repulsion integrals
- Electron density on grids
- Electrostatic potential
- Moments (dipole, quadrupole, etc.)
- Basis set parsing
- General contractions

**Sources**: GitHub repository, QCDevs project, ReadTheDocs

## Key Strengths

### Modular Design:
- Separate evaluation and integration
- Composable components
- Extensible architecture
- Clean API

### Pure Python:
- NumPy/SciPy based
- No compilation required
- Cross-platform
- Easy installation

### QCDevs Ecosystem:
- Part of larger project
- Integration with other tools
- Maintained community
- Educational focus

### Comprehensive Integrals:
- All standard molecular integrals
- Grid-based quantities
- Properties and moments
- Flexible basis handling

## Inputs & Outputs
- **Input formats**:
  - IOData molecule objects
  - Basis set dictionaries
  - NumPy arrays for coordinates
  
- **Output data types**:
  - Integral arrays
  - Density values
  - Potential values
  - Moment tensors

## Interfaces & Ecosystem
- **QCDevs tools**:
  - IOData (file I/O)
  - Grid (numerical grids)
  - QCSchema compatibility
  
- **External**:
  - NumPy/SciPy
  - Standard basis set formats
  - Interoperable with other codes

## Advanced Features

### Integral Evaluation:
- Cartesian and spherical GTOs
- Contracted basis functions
- General contraction support
- Efficient algorithms

### Property Calculations:
- Multipole moments
- Electrostatic potential
- Electron density
- Gradient quantities

### Basis Set Handling:
- Standard library formats
- Custom basis sets
- Segmented and general contractions

## Performance Characteristics
- **Speed**: NumPy vectorized
- **Accuracy**: Numerical precision
- **System size**: Small to medium
- **Memory**: NumPy arrays
- **Parallelization**: NumPy threading

## Computational Cost
- **Integrals**: Standard N^4 for ERIs
- **Grids**: Linear in grid points
- **Typical**: Development/research scale

## Limitations & Known Constraints
- **Not full QM code**: Library, not standalone
- **Large systems**: Not optimized for production
- **Speed**: Python overhead
- **Features**: Core integrals focus

## Comparison with Other Codes
- **vs libint**: GBasis Python, libint C++
- **vs PySCF integrals**: Different scope
- **Unique strength**: Pure Python, educational, modular

## Application Areas

### Code Development:
- New QM method implementation
- Testing integral routines
- Algorithm development
- Prototyping

### Education:
- Teaching molecular integrals
- Understanding QM code internals
- Research training

### Analysis:
- Density evaluation
- Property calculations
- Wavefunction analysis

### Integration:
- Building custom QM codes
- Interfacing with other tools
- Workflow development

## Best Practices

### Getting Started:
- pip install gbasis
- Use IOData for molecules
- Follow documentation examples

### Integral Computation:
- Select appropriate functions
- Understand basis conventions
- Validate against references

## Community and Support
- Open source LGPL v3
- QCDevs community
- GitHub development
- ReadTheDocs documentation
- PyPI distribution

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/theochem/gbasis
2. Documentation: https://gbasis.readthedocs.io/
3. PyPI: https://pypi.org/project/gbasis/
4. QCDevs project

**Confidence**: VERIFIED - Active development, published

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, LGPL v3)
- Package: PyPI
- Documentation: ReadTheDocs
- Community: QCDevs
- Specialty: Gaussian integral library, modular design
