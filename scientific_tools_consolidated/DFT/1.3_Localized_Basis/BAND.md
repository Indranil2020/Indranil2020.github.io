# BAND (Amsterdam Modeling Suite)

## Official Resources
- Homepage: https://www.scm.com/product/band/
- Documentation: https://www.scm.com/doc/BAND/
- Tutorials: https://www.scm.com/doc/Tutorials/BAND/
- Developer: Software for Chemistry & Materials (SCM)
- License: Commercial (Academic pricing available)

## Overview
BAND is the periodic Density Functional Theory (DFT) code within the Amsterdam Modeling Suite (AMS). Unlike most periodic codes that use plane waves (like VASP or QE), BAND utilizes atom-centered numerical orbitals (STOs/NAOs). This basis set allows for an accurate treatment of both core and valence electrons and makes the code particularly efficient for low-dimensional systems (1D polymers, 2D slabs) and empty space.

**Scientific domain**: Periodic systems, surface science, catalysis, nanotubes, polymers  
**Target user community**: Chemists and materials scientists studying periodic systems where chemical insight (orbitals, bonds) is paramount

## Theoretical Methods
- Density Functional Theory (DFT)
- Slater-Type Orbitals (STOs) and Numerical Atomic Orbitals (NAOs)
- LDA, GGA, meta-GGA exchange-correlation functionals
- Hybrid functionals (B3LYP, PBE0, HSE)
- Dispersion corrections (DFT-D3, DFT-D4)
- Scalar relativistic ZORA
- Spin-orbit coupling
- COSMO solvation for periodic surfaces
- DFT+U for correlated systems

## Capabilities (CRITICAL)
- Ground-state electronic structure
- 1D periodic (nanowires, polymers)
- 2D periodic (surfaces, slabs)
- 3D periodic (bulk crystals)
- Band structure and DOS
- COOP/COHP bonding analysis
- Geometry optimization
- Transition state search (NEB)
- Phonon calculations
- EELS spectra
- STM image simulation
- Work function calculations
- Forces and stress tensors

**Sources**: SCM official documentation, AMS Suite manuals

## Key Strengths

### Slater-Type Orbitals:
- Correct electron cusp behavior
- Compact basis representation
- Chemical orbital interpretation
- Efficient for open structures
- Better asymptotic decay

### True Low-Dimensional Periodicity:
- Native 1D/2D periodicity
- No artificial vacuum padding
- Correct electrostatics
- Efficient for surfaces
- Polymer chain calculations

### Chemical Bonding Analysis:
- COOP (Crystal Orbital Overlap Population)
- COHP (Crystal Orbital Hamilton Population)
- Mulliken analysis
- Hirshfeld charges
- Band decomposition

### Relativistic Treatment:
- ZORA (scalar relativistic)
- Spin-orbit coupling
- Heavy element support
- Accurate for actinides
- Core electron treatment

## Inputs & Outputs
- **Input formats**:
  - AMS GUI structure builder
  - CIF files
  - XYZ with lattice vectors
  - POSCAR (via conversion)
  
- **Output data types**:
  - Total energies
  - Band structures
  - DOS/PDOS
  - COOP/COHP data
  - Optimized structures
  - Phonon spectra
  - STM images

## Interfaces & Ecosystem
- **AMS Integration**:
  - Unified driver for all SCM codes
  - Seamless ADF-BAND coupling
  - ReaxFF interface
  - DFTB interface
  
- **Python scripting**:
  - PLAMS (Python Library for AMS)
  - Workflow automation
  - Batch processing
  - Custom analysis
  
- **Visualization**:
  - AMS-GUI (native)
  - ADFView for orbitals
  - Band structure plotter
  - DOS visualization

## Advanced Features

### Fragment Analysis:
- Periodic fragment calculations
- Adsorbate-surface decomposition
- Energy decomposition analysis
- Bonding contributions

### Surface Science:
- Work function calculations
- Adsorbate binding energies
- Surface reconstruction
- Step edge modeling

### Hybrid Functionals:
- HSE06 for band gaps
- PBE0 calculations
- Range-separated hybrids
- Accurate gap prediction

### Phonon Calculations:
- Finite differences
- IR intensities
- Thermodynamic properties
- Free energy calculations

### QM/MM Coupling:
- Periodic QM/MM
- Embedding in force fields
- Multi-scale modeling

## Performance Characteristics
- **Speed**: Efficient for open structures
- **Accuracy**: High with appropriate basis
- **System size**: Hundreds of atoms typical
- **Memory**: Basis-dependent
- **Parallelization**: Hybrid MPI/OpenMP

## Computational Cost
- **Basis scaling**: TZP costs ~3x DZ
- **Meta-GGA**: 2-4x GGA cost
- **SOC**: 4-8x scalar relativistic
- **Hybrids**: 10-50x GGA cost
- **Typical**: Workstation to cluster

## Limitations & Known Constraints
- **Commercial license**: Required for use
- **Basis convergence**: Large basis needed for high accuracy
- **Linear dependence**: Issue with dense large basis sets
- **Metallic systems**: Can be challenging
- **Hybrid cost**: Expensive for large systems

## Comparison with Other Codes
- **vs VASP/QE**: BAND uses STOs, plane-wave codes use PWs; BAND better for open structures
- **vs CRYSTAL**: Both localized basis; BAND STOs, CRYSTAL GTOs
- **vs SIESTA**: Both NAO-based; different commercial/open model
- **Unique strength**: STOs, true 1D/2D periodicity, COOP/COHP analysis, AMS integration

## Application Areas

### Surface Science:
- Catalytic surfaces
- Adsorbate binding
- Surface reconstruction
- Work function engineering

### Nanomaterials:
- Carbon nanotubes
- Graphene nanoribbons
- 2D materials
- Nanowires

### Polymers:
- Conjugated polymers
- Polymer electronics
- Band engineering
- Conducting polymers

### Heavy Elements:
- Lanthanide/actinide compounds
- Relativistic effects
- f-electron systems
- Nuclear materials

## Best Practices

### Basis Set Selection:
- Start with DZ/TZP for geometry
- Use QZ4P for final properties
- Check linear dependency
- Test basis convergence

### SCF Convergence:
- Use NumericalAccuracy setting
- DIIS for difficult cases
- Level shifting if needed
- Check k-point convergence

### Relativistic Calculations:
- ZORA scalar as default
- SOC only when needed
- Core treatment consistent
- Heavy elements require testing

### Performance:
- Pure MPI for multi-node
- OpenMP for single node
- Balance k-points and cores
- Monitor memory usage

## Community and Support
- Professional SCM support
- Extensive documentation
- Video tutorials
- Training workshops
- Active user forum

## Verification & Sources
**Primary sources**:
1. Official website: https://www.scm.com/product/band/
2. Documentation: https://www.scm.com/doc/BAND/
3. G. te Velde, E.J. Baerends et al., J. Comput. Chem. publications

**Secondary sources**:
1. SCM tutorials
2. Published applications
3. Benchmark studies

**Confidence**: CONFIRMED - Commercial product, established code

**Verification status**: ✅ VERIFIED
- Source code: Commercial (SCM)
- Documentation: Extensive
- Support: Professional
- Active development: Regular releases
- Specialty: STO periodic DFT, true 1D/2D, COOP/COHP analysis

