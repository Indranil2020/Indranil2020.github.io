# HORTON (Helpful Open-source Research TOol for N-fermion systems)

## Official Resources
- Homepage: https://theochem.github.io/horton/
- Documentation: https://horton.readthedocs.io/
- Source Repository: https://github.com/theochem/horton
- License: GNU General Public License v3.0

## Overview
HORTON is a modular quantum chemistry program written primarily in Python, designed for electronic structure calculations, method prototyping, and educational purposes. It emphasizes code readability, extensibility, and user-friendliness over raw computational performance. HORTON 3.x has evolved into a suite of independent modules (GBasis, Grid, IOData) that work together to provide a flexible framework for quantum chemistry workflows.

**Scientific domain**: Molecules, quantum chemistry, method development, education  
**Target user community**: Researchers developing new methods, educators teaching quantum chemistry, those needing interpretable and extensible code

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Density Functional Theory (DFT)
- Gaussian Type Orbitals (GTOs)
- Exchange-correlation via LibXC
- Post-HF methods (MP2, limited)
- Orbital optimization
- Density matrix methods
- Conceptual DFT descriptors

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Hartree-Fock calculations
- DFT with multiple functionals
- Molecular integrals (via GBasis)
- Numerical integration grids (via Grid)
- File I/O for multiple formats (via IOData)
- Orbital localization
- Population analysis
- Conceptual DFT reactivity descriptors
- Wavefunction analysis
- Density-based properties

**Sources**: GitHub repository, QCDevs community, publications

## Key Strengths

### Modular Architecture:
- GBasis: Gaussian integral evaluation
- Grid: Numerical integration grids
- IOData: File format parsing
- Independent, reusable components
- Mix-and-match functionality

### Python-Native Design:
- Pure Python with NumPy/SciPy
- Readable, educational code
- Interactive development
- Jupyter notebook compatible
- Easy to extend and modify

### Conceptual DFT:
- Reactivity descriptors
- Fukui functions
- Hardness and softness
- Electrophilicity indices
- Dual descriptors

### Educational Focus:
- Transparent algorithms
- Well-documented code
- Teaching-oriented design
- Prototyping-friendly
- Method development platform

## Inputs & Outputs
- **Input formats**:
  - XYZ coordinates
  - Gaussian fchk files
  - Molden files
  - WFN/WFX files
  - Multiple QC output formats (via IOData)
  
- **Output data types**:
  - Total energies
  - Orbital data
  - Density matrices
  - Molecular properties
  - Analysis results

## Interfaces & Ecosystem
- **QCDevs project**:
  - GBasis for integrals
  - Grid for numerical grids
  - IOData for file I/O
  - ChemTools for analysis
  
- **External integration**:
  - LibXC for functionals
  - NumPy/SciPy ecosystem
  - Matplotlib visualization
  - Jupyter notebooks

## Advanced Features

### GBasis Module:
- One-electron integrals
- Two-electron integrals
- Molecular orbital evaluation
- Density evaluation on grids
- Property calculations

### Grid Module:
- Becke partitioning
- Atom-centered grids
- Lebedev angular grids
- Radial grid schemes
- Integration accuracy control

### IOData Module:
- Read/write 40+ file formats
- Quantum chemistry packages
- Molecular dynamics formats
- Plane-wave DFT outputs
- Format conversion utilities

### Conceptual DFT:
- Local reactivity indices
- Global descriptors
- Condensed-to-atoms
- Orbital-based analysis

## Performance Characteristics
- **Speed**: Educational, not production-optimized
- **Accuracy**: Standard quantum chemistry
- **System size**: Small to medium molecules
- **Memory**: Python/NumPy requirements
- **Parallelization**: Limited (NumPy threading)

## Computational Cost
- **Focus**: Understanding over speed
- **Typical**: Seconds to minutes for small molecules
- **Purpose**: Prototyping and education
- **Scaling**: Standard Gaussian basis scaling

## Limitations & Known Constraints
- **Production use**: Not for production calculations
- **System size**: Small molecules only
- **Speed**: Slower than compiled codes
- **Periodicity**: Molecular focus
- **Post-HF**: Limited beyond MP2
- **Active development**: HORTON 3 modular rewrite ongoing

## Comparison with Other Codes
- **vs PySCF**: HORTON more modular, PySCF more complete
- **vs PyDFT**: Similar educational goals, different scope
- **vs Gaussian**: HORTON open, educational; Gaussian production
- **Unique strength**: Modular design, conceptual DFT, educational focus, IOData ecosystem

## Application Areas

### Education:
- Teaching quantum chemistry
- Algorithm understanding
- Graduate courses
- Self-study
- Computational chemistry workshops

### Method Development:
- New functional testing
- Algorithm prototyping
- Method validation
- Research exploration
- Conceptual DFT research

### Wavefunction Analysis:
- Population analysis
- Orbital localization
- Density partitioning
- Bonding characterization
- Reactivity prediction

### File Processing:
- Format conversion (IOData)
- Data extraction
- Workflow integration
- Multi-code projects

## Best Practices

### Getting Started:
- Install individual modules (gbasis, grid, iodata)
- Follow documentation tutorials
- Start with small test systems
- Use Jupyter for exploration

### Development:
- Leverage modular components
- Write readable code
- Document thoroughly
- Test against references

### Integration:
- Use IOData for file I/O
- Combine with other tools
- Script complex workflows
- Validate results

## Community and Support
- Open source GPL v3
- QCDevs community
- GitHub repositories
- ReadTheDocs documentation
- Academic publications
- Growing contributor base

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/theochem/horton
2. QCDevs: https://qcdevs.org/
3. T. Verstraelen et al., J. Chem. Theory Comput. publications

**Secondary sources**:
1. GBasis, Grid, IOData papers
2. Conceptual DFT literature
3. Python quantum chemistry surveys

**Confidence**: VERIFIED - Active development, published methodology

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Academic use: Published applications
- Documentation: ReadTheDocs
- Active development: Modular rewrite ongoing
- Specialty: Modular Python QC, conceptual DFT, education, IOData ecosystem
