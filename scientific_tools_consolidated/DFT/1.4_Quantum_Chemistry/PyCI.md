# PyCI

## Official Resources
- Homepage: https://github.com/theochem/pyci
- Documentation: https://pyci.readthedocs.io/
- Source Repository: https://github.com/theochem/pyci
- License: GNU General Public License v3.0

## Overview
PyCI is a Python library for configuration interaction (CI) calculations. Part of the TheoChem ecosystem, it provides efficient CI implementations using determinant-based algorithms with optimized Slater-Condon rules. Integrates with ModelHamiltonian and FanPy for comprehensive wavefunction studies.

**Scientific domain**: Configuration interaction, electron correlation, FCI  
**Target user community**: Researchers performing CI calculations and wavefunction method development

## Theoretical Methods
- Full Configuration Interaction (FCI)
- Selected Configuration Interaction (sCI)
- CISD (singles and doubles)
- CISDT (singles, doubles, triples)
- Active space CI
- Heat-bath CI (HCI)
- Spin-adapted implementations
- Sparse Hamiltonian methods

## Capabilities (CRITICAL)
- Full CI for small systems
- Selected CI for larger systems
- Determinant-based algorithms
- Efficient Slater-Condon rules
- Sparse matrix techniques
- ModelHamiltonian integration
- FanPy wavefunction interface
- Custom Hamiltonian support
- Ground and excited states
- Spin eigenfunctions

## Key Strengths

### CI Methods:
- Multiple truncation levels
- Selection schemes (HCI)
- Variational treatment
- Size-consistency corrections

### Efficiency:
- Optimized Slater-Condon
- Sparse Hamiltonian storage
- Davidson diagonalization
- Memory-efficient schemes

### Integration:
- TheoChem ecosystem
- Standard integral formats
- FanPy wavefunctions
- Model Hamiltonians

### Research Flexibility:
- Custom CI spaces
- User-defined selections
- Algorithm testing
- Method development

## Inputs & Outputs
- **Input formats**:
  - Python API
  - FCIDUMP integrals
  - ModelHamiltonian objects
  
- **Output data types**:
  - CI energies
  - CI vectors
  - Natural orbitals
  - Excited states

## Interfaces & Ecosystem
- **TheoChem tools**: FanPy, ModelHamiltonian
- **Integral sources**: PySCF interface
- **NumPy/SciPy**: Sparse linear algebra
- **Diagonalization**: Davidson, Lanczos

## Advanced Features

### Selected CI:
- Heat-bath selection
- Importance truncation
- Perturbative corrections
- Convergence extrapolation

### Spin Adaptation:
- Sz eigenfunctions
- S² eigenfunctions
- Genealogical coupling
- CSF basis

### Sparse Methods:
- On-the-fly generation
- Compressed storage
- Efficient contractions
- Large CI spaces

## Performance Characteristics
- **Speed**: Efficient determinant handling
- **Accuracy**: Variational CI accuracy
- **System size**: Small to medium (FCI), larger (sCI)
- **Memory**: Sparse storage efficient
- **Parallelization**: Potential for parallel

## Computational Cost
- **FCI**: Factorial scaling (small systems)
- **CISD**: O(N^6) manageable
- **Selected CI**: System dependent
- **Typical**: Seconds to hours

## Limitations & Known Constraints
- **FCI scaling**: Exponential, limits system size
- **Size consistency**: Not inherent in truncated CI
- **Production**: Research focus
- **Documentation**: Academic oriented

## Comparison with Other Codes
- **vs MOLPRO CI**: PyCI more flexible, MOLPRO faster
- **vs PySCF CI**: Similar, different ecosystems
- **vs Arrow/DICE**: Different sCI algorithms
- **Unique strength**: TheoChem integration, flexibility

## Application Areas

### Benchmarks:
- FCI reference energies
- Correlation energy studies
- Method validation

### Strong Correlation:
- Multi-reference systems
- Bond breaking
- Transition metals

### Method Development:
- CI algorithm testing
- Selection scheme development
- Hybrid methods

## Best Practices

### Calculation Setup:
- Appropriate active space
- Selection thresholds
- Convergence criteria
- State averaging for excited states

### Efficiency:
- Start with smaller CI
- Increase systematically
- Use selected CI for larger systems
- Monitor memory usage

## Community and Support
- Open-source GPL v3
- TheoChem group (McMaster University)
- Academic publications
- GitHub for issues
- Growing user base

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/theochem/pyci
2. Ayers group publications
3. CI methodology papers
4. TheoChem documentation

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Documentation: ReadTheDocs
- Academic group: TheoChem
- Active development: Yes
