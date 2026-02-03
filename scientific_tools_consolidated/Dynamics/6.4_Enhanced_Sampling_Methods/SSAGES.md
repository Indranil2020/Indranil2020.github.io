# SSAGES

## Official Resources
- Homepage: https://ssagesproject.github.io/
- Documentation: https://ssages.readthedocs.io/
- Source Repository: https://github.com/SSAGESproject/SSAGES
- License: GPL-3.0

## Overview
SSAGES (Software Suite for Advanced General Ensemble Simulations) is a free, open-source software package for performing advanced sampling simulations. It provides a unified interface to multiple enhanced sampling methods and integrates with popular MD engines.

**Scientific domain**: Enhanced sampling, free energy calculations, rare events  
**Target user community**: Researchers studying rare events and free energy landscapes

## Theoretical Methods
- Metadynamics
- Adaptive biasing force (ABF)
- Umbrella sampling
- Forward flux sampling
- String method
- Basis function sampling

## Capabilities (CRITICAL)
- Multiple enhanced sampling methods
- Multiple MD engine support
- Collective variable library
- Free energy calculations
- Rare event sampling
- C++ implementation

## Key Strengths

### Method Variety:
- Many enhanced sampling methods
- Unified interface
- Easy method switching
- Extensible

### MD Engine Support:
- LAMMPS
- GROMACS
- OpenMD
- Hoomd-blue
- QBox

## Inputs & Outputs
- **Input formats**:
  - JSON configuration
  - MD engine inputs
  
- **Output data types**:
  - Free energy profiles
  - Collective variable trajectories
  - Bias potentials

## Interfaces & Ecosystem
- **LAMMPS**: Integration
- **GROMACS**: Integration
- **HOOMD-blue**: Integration
- **PySAGES**: Python version

## Advanced Features
- **ABF**: Adaptive biasing force
- **Metadynamics**: History-dependent bias
- **String method**: Reaction pathways
- **FFS**: Forward flux sampling
- **Custom CVs**: User-defined variables

## Performance Characteristics
- C++ implementation
- Efficient CV calculation
- Good parallel scaling
- Low overhead

## Computational Cost
- Overhead depends on method
- ABF/metadynamics moderate
- String method more expensive
- Overall: Efficient

## Best Practices
- Choose appropriate method
- Validate CV choice
- Check convergence
- Use sufficient sampling

## Limitations & Known Constraints
- C++ complexity
- Setup can be involved
- Documentation varies
- PySAGES easier to use

## Application Areas
- Protein folding
- Chemical reactions
- Phase transitions
- Nucleation
- Conformational changes

## Comparison with Other Codes
- **vs PLUMED**: SSAGES C++ standalone, PLUMED plugin architecture
- **vs PySAGES**: SSAGES C++, PySAGES Python/GPU
- **vs Colvars**: SSAGES more methods, Colvars more CV types
- **Unique strength**: Unified interface to many methods, multiple MD engine support

## Community and Support
- Active development
- GitHub issues
- Documentation
- PySAGES Python alternative

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/SSAGESproject/SSAGES
2. H. Sidky et al., J. Chem. Phys. 148, 044104 (2018)

**Secondary sources**:
1. SSAGES tutorials
2. PySAGES documentation
3. Enhanced sampling publications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-3.0)
- Academic citations: >200
- Active development
