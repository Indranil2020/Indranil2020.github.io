# FanPy

## Official Resources
- Homepage: https://github.com/theochem/fanpy
- Documentation: https://fanpy.readthedocs.io/
- Source Repository: https://github.com/theochem/fanpy
- License: GNU General Public License v3.0

## Overview
FanPy (Flexible Ansatz for N-electron Wavefunctions in Python) is a Python library for ab initio quantum chemistry calculations using flexible wavefunction ansätze. It enables development and application of novel correlated wavefunction methods, particularly geminal-based approaches. Part of the TheoChem ecosystem alongside ModelHamiltonian and PyCI.

**Scientific domain**: Correlated wavefunction methods, geminal theory, novel ansätze  
**Target user community**: Researchers developing and applying novel wavefunction methods for strong correlation

## Theoretical Methods
- Antisymmetric Product of Geminals (APG)
- Antisymmetrized Product of Strongly Orthogonal Geminals (APSG)
- Antisymmetric Product of Interacting Geminals (APIG)
- Generalized valence bond (GVB) separations
- Configuration interaction
- Coupled cluster variants
- Custom wavefunction ansätze
- Parametrized wavefunction optimization

## Capabilities (CRITICAL)
- Flexible wavefunction ansätze
- Geminal-based methods for strong correlation
- Variational optimization of wavefunction parameters
- ProjectQ integration for custom operators
- ModelHamiltonian integration for model systems
- Linear and non-linear parameter optimization
- Ground and excited state calculations
- Research and development platform
- Extensive customization options

## Key Strengths

### Flexibility:
- Custom ansätze design
- Novel wavefunction forms
- Research-oriented architecture
- Easy method prototyping
- Modular components

### Geminal Methods:
- APG wavefunctions
- APSG for strong correlation
- APIG with geminal interactions
- Seniority-zero methods
- Compact correlation descriptions

### Strong Correlation:
- Bond breaking scenarios
- Transition metal complexes
- Biradicals and polyradicals
- Multi-reference character

### Integration:
- ModelHamiltonian for model systems
- PyCI for CI wavefunctions
- TheoChem ecosystem
- Standard integral formats

## Inputs & Outputs
- **Input formats**:
  - Python API
  - Integral files (fcidump format)
  - ModelHamiltonian objects
  
- **Output data types**:
  - Optimized energies
  - Wavefunction parameters
  - Geminal coefficients
  - Occupation numbers

## Interfaces & Ecosystem
- **TheoChem tools**: ModelHamiltonian, PyCI
- **Integral packages**: PySCF, horton
- **NumPy/SciPy**: Numerical operations
- **Optimization**: Various optimizers supported

## Advanced Features

### Geminal Optimization:
- Orbital optimization
- Geminal coefficient optimization
- Coupled optimization schemes
- Symmetry constraints

### Custom Ansätze:
- User-defined wavefunction forms
- Operator algebra
- Flexible parameterization
- Research development

### Model Systems:
- Hubbard models
- PPP models
- Custom Hamiltonians
- Lattice systems

## Performance Characteristics
- **Speed**: Python implementation, moderate
- **Accuracy**: Depends on ansatz choice
- **System size**: Small to medium molecules
- **Memory**: Scales with ansatz complexity
- **Parallelization**: NumPy/SciPy threading

## Computational Cost
- **APG**: Polynomial scaling
- **APSG**: Reduced active space
- **Optimization**: Depends on parameters
- **Typical**: Minutes to hours for small systems

## Limitations & Known Constraints
- **System size**: Best for small molecules
- **Production**: Research-focused
- **Documentation**: Academic oriented
- **Community**: Research group centered
- **Learning curve**: Requires QC background

## Comparison with Other Codes
- **vs PySCF**: FanPy specialized in geminals
- **vs Molpro GVB**: FanPy more flexible
- **vs PyCI**: Different focus (geminals vs CI)
- **Unique strength**: Novel geminal methods, flexibility

## Application Areas

### Strong Correlation:
- Bond dissociation
- Transition metal complexes
- Open-shell systems
- Multi-reference problems

### Method Development:
- New ansätze testing
- Algorithm development
- Proof of concept
- Comparison studies

### Model Systems:
- Hubbard model studies
- Lattice models
- Strongly correlated models
- Benchmark calculations

## Best Practices

### Ansatz Selection:
- Match ansatz to problem
- Start simple, increase complexity
- Validate with known results
- Monitor convergence

### Optimization:
- Choose appropriate optimizer
- Reasonable initial guess
- Convergence thresholds
- Multiple restarts if needed

## Community and Support
- Open-source GPL v3
- TheoChem group (McMaster University)
- Academic publications
- GitHub issues for support
- Growing community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/theochem/fanpy
2. Ayers group publications
3. Geminal method papers
4. TheoChem ecosystem documentation

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Documentation: ReadTheDocs
- Academic group: TheoChem, McMaster
- Active development: Yes
