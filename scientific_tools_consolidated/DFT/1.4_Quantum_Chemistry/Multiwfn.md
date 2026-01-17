# Multiwfn

## Official Resources
- Homepage: http://sobereva.com/multiwfn/
- Documentation: http://sobereva.com/multiwfn/Multiwfn_manual.html
- Download: http://sobereva.com/multiwfn/download.html
- License: Free for academic use

## Overview
Multiwfn is a comprehensive, extremely powerful electron wavefunction analysis toolbox. It can perform a wide variety of wavefunction analyses based on output from almost all major quantum chemistry packages including Gaussian, ORCA, GAMESS, NWChem, Molpro, and many others.

**Scientific domain**: Wavefunction analysis, bonding analysis, property calculations  
**Target user community**: Computational chemists needing detailed wavefunction and property analysis

## Theoretical Methods (Analysis)
- Quantum Theory of Atoms in Molecules (QTAIM)
- Natural Bond Orbital analysis (NBO-like)
- Electron Localization Function (ELF)
- Localized Orbital Locator (LOL)
- Reduced Density Gradient (NCI)
- Electron Density Difference
- Orbital composition analysis
- Population analysis methods
- Excited state analysis
- Aromaticity indices

## Capabilities (CRITICAL)
- QTAIM topology analysis
- Bond critical point analysis
- Basin integration
- ELF/LOL visualization
- Molecular surface analysis
- ESP mapping
- Hirshfeld/CM5 charges
- Mayer bond orders
- Natural atomic orbitals
- TDDFT analysis
- Hole-electron analysis
- Absorption/emission spectra
- Fukui functions
- 3D grid calculations

## Key Strengths

### Universal Input:
- Gaussian output files
- ORCA output/molden
- GAMESS output
- NWChem output
- Molpro output
- PSI4 output
- Q-Chem output
- Generic wfn/wfx/fchk

### Analysis Breadth:
- Hundreds of analysis functions
- Property visualization
- Quantitative metrics
- Publication-ready output

### User Interface:
- Interactive menus
- Batch processing
- Scripting capability
- Graphical output

### Documentation:
- Comprehensive manual (800+ pages)
- Tutorials
- Active support
- Regular updates

## Inputs & Outputs
- **Input formats**:
  - fchk (Gaussian)
  - molden
  - wfn/wfx
  - cube files
  - NBO output
  - Various program outputs
  
- **Output data types**:
  - Critical point data
  - Basin properties
  - Population charges
  - Spectral data
  - 3D visualization files
  - Publication tables

## Interfaces & Ecosystem
- **Visualization**: VMD, GaussView, Chemcraft
- **QC codes**: Gaussian, ORCA, GAMESS, NWChem, Molpro, etc.
- **File formats**: Cube, molden, wfn/wfx, xyz
- **Post-processing**: Scripting, batch mode

## Advanced Features

### QTAIM Analysis:
- Automated CP search
- Basin integration
- IQA energy decomposition
- Source function
- Delocalization indices

### Bonding Analysis:
- ELF/LOL
- NCI analysis
- Bond order calculations
- Orbital contributions
- Fragment analyses

### Property Mapping:
- Electrostatic potential
- Fukui functions
- ALIE
- Electron density
- Custom properties

### Spectroscopy:
- UV-Vis spectra
- Emission spectra
- Vibrational analysis
- NMR predictions
- Circular dichroism

## Performance Characteristics
- **Speed**: Efficient for analysis
- **Accuracy**: High-precision integration
- **System size**: Thousands of atoms
- **Memory**: Manageable
- **Platform**: Windows/Linux/macOS

## Computational Cost
- **QTAIM**: Moderate (integration)
- **Grid analysis**: Grid-size dependent
- **Batch mode**: Scriptable
- **Large systems**: Efficient algorithms
- **Typical**: Minutes for molecular analysis

## Limitations & Known Constraints
- **Analysis only**: No electronic structure calculation
- **Input dependent**: Quality depends on source calculation
- **Learning curve**: Many options require understanding
- **Visualization**: External tools for 3D
- **Closed source**: Free but not open

## Comparison with Other Codes
- **vs AIMALL**: Both QTAIM; Multiwfn broader, AIMALL specialized
- **vs NBO**: Different focus; Multiwfn more diverse
- **vs Critic2**: Multiwfn more features, Critic2 periodic focus
- **vs Postprocessors**: Most comprehensive analysis tool
- **Unique strength**: Breadth of methods, universal input

## Application Areas

### Chemical Bonding:
- Bond characterization
- Weak interactions (NCI)
- Aromaticity
- Hyperconjugation

### Reactivity:
- Fukui function analysis
- Electrophilicity
- Local reactivity descriptors
- Reaction mechanisms

### Spectroscopy:
- Band assignment
- Transition analysis
- Hole-electron distributions
- Optical properties

### Materials:
- Surface properties
- Charge distribution
- Intermolecular interactions
- Crystal analysis

## Best Practices

### Input Preparation:
- Appropriate basis sets
- Sufficient grid density
- Proper wavefunction quality
- Complete output saving

### Analysis Selection:
- Match method to question
- Understand limitations
- Validate with multiple methods
- Check convergence

## Community and Support
- Free for academic use
- Active development (T. Lu)
- Extensive manual
- Email support
- Publication citations (5000+)

## Verification & Sources
**Primary sources**:
1. Official website: http://sobereva.com/multiwfn/
2. Lu, Chen, J. Comput. Chem. 33, 580-592 (2012)
3. Manual: 800+ page documentation
4. Regular updates and releases

**Confidence**: CONFIRMED
- Software: Distributed free for academics
- Documentation: Comprehensive
- Citations: >5000 in literature
- Active development: Yes
- Community: Large, active user base
