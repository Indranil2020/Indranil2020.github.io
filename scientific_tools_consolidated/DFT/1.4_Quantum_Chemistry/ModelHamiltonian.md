# ModelHamiltonian

## Official Resources
- Homepage: https://github.com/theochem/ModelHamiltonian
- Documentation: https://modelhamiltonian.readthedocs.io/
- Source Repository: https://github.com/theochem/ModelHamiltonian
- License: GNU General Public License v3.0

## Overview
ModelHamiltonian is a Python library that facilitates the application of quantum chemistry methods to model Hamiltonians by translating them into standard 0-, 1-, and 2-electron integrals. It bridges model system studies with ab initio quantum chemistry codes, enabling use of sophisticated wavefunction methods on lattice models.

**Scientific domain**: Model Hamiltonians, condensed matter physics, quantum chemistry bridges  
**Target user community**: Researchers studying model systems with quantum chemistry methods

## Theoretical Methods
- Hubbard model (various geometries)
- Pariser-Parr-Pople (PPP) model
- Heisenberg model mapping
- Anderson impurity model
- Extended Hubbard models
- Custom model Hamiltonians
- Lattice geometries (1D, 2D, 3D)

## Capabilities (CRITICAL)
- Model to integral conversion
- Standard 0/1/2-electron format
- FCIDUMP output format
- FanPy wavefunction integration
- PyCI CI integration
- Custom Hamiltonians definition
- Parameter specification (t, U, V, J)
- Lattice model generation
- Periodic and open boundaries
- Various lattice geometries

## Key Strengths

### Model Translation:
- Standard integral format
- Interoperability with any QC code
- Easy FCIDUMP generation
- Flexible parameter specification

### Lattice Models:
- 1D chains
- 2D square/triangular/honeycomb
- 3D cubic systems
- Custom topologies

### Physical Models:
- Hubbard for correlation
- PPP for π-conjugated systems
- Extended models (V, t')
- Periodic Anderson

### Ecosystem:
- FanPy for geminal methods
- PyCI for CI calculations
- Standard QC code interfaces
- TheoChem tool integration

## Inputs & Outputs
- **Input formats**:
  - Python API
  - Lattice specifications
  - Parameter dictionaries
  
- **Output data types**:
  - FCIDUMP files
  - Integral arrays (0/1/2 electron)
  - Hamiltonian matrices
  - Lattice visualizations

## Interfaces & Ecosystem
- **TheoChem tools**: FanPy, PyCI
- **Output formats**: FCIDUMP standard
- **NumPy/SciPy**: Array handling
- **Visualization**: Lattice plotting

## Advanced Features

### Hubbard Models:
- On-site interaction (U)
- Hopping (t)
- Next-nearest neighbor (t')
- Extended interactions (V)

### PPP Model:
- Ohno potential
- Mataga-Nishimoto
- Custom screening
- π-electron systems

### Custom Hamiltonians:
- User-defined terms
- Arbitrary operators
- Parameter sweeps
- Model development

### Geometry Support:
- Chain, ring, ladder
- Square, triangular, honeycomb
- Bethe lattice
- Custom connectivity

## Performance Characteristics
- **Speed**: Fast generation
- **Accuracy**: Exact model representation
- **System size**: Limited by subsequent QC
- **Memory**: Scales with lattice size
- **Output**: Immediate generation

## Computational Cost
- **Generation**: Milliseconds to seconds
- **Bottleneck**: Subsequent QC calculation
- **Large lattices**: Memory for integrals
- **Typical**: Fast preprocessing step

## Limitations & Known Constraints
- **Model focus**: Not for real materials
- **Ab initio**: No molecular calculations
- **Physical insight**: Requires model understanding
- **Size**: QC method limits system size

## Comparison with Other Codes
- **vs Direct Hubbard codes**: More flexible output
- **vs DMRG codes**: Different focus (bridging)
- **vs PySCF model**: Similar, different ecosystem
- **Unique strength**: QC method interoperability

## Application Areas

### Strongly Correlated Systems:
- Mott insulators
- Antiferromagnetism
- Superconductivity pairing
- Quantum phase transitions

### Conjugated Systems:
- Polyenes
- Graphene fragments
- Organic semiconductors
- π-electron models

### Method Development:
- Algorithm testing
- Method validation
- New ansätze testing
- Comparison studies

### Education:
- Teaching correlation
- Model system exploration
- Visualization of physics
- Student projects

## Best Practices

### Model Selection:
- Match physics to model
- Appropriate parameters
- Validate with literature
- Size convergence

### Parameter Choice:
- Physical values from literature
- Sensitivity analysis
- Multiple U/t ratios
- Half-filling studies

## Community and Support
- Open-source GPL v3
- TheoChem group (McMaster University)
- Academic publications
- GitHub for contributions
- Documentation and examples

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/theochem/ModelHamiltonian
2. Ayers group publications
3. Hubbard/PPP model literature
4. TheoChem ecosystem documentation

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Documentation: ReadTheDocs
- Academic group: TheoChem
- Active development: Yes
