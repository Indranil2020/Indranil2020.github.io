# AMPAC

## Official Resources
- Homepage: https://openmopac.net/ (MOPAC as successor)
- Documentation: Historic manuals
- Note: Austin-Method Package
- License: Academic/Commercial (historic)

## Overview
AMPAC (Austin-Method Package) is a historic semi-empirical quantum chemistry program package implementing AM1, PM3, and related methods. Developed with significant contributions from Michael Dewar's group at the University of Texas Austin, it represented a major advance in semi-empirical methodology for organic chemistry applications.

**Scientific domain**: Semi-empirical quantum chemistry  
**Target user community**: Historic; mostly superseded by MOPAC, but methods remain foundational

## Theoretical Methods
- AM1 (Austin Model 1)
- PM3 (Parametric Method 3)
- MNDO (Modified Neglect of Diatomic Overlap)
- SAM1 (Semi-Ab initio Model 1)
- Geometry optimization
- Transition state searches
- Frequency calculations
- Reaction path following

## Capabilities (CRITICAL)
- Fast semi-empirical energies
- Full geometry optimization
- Transition state location
- Vibrational frequencies
- Thermodynamic properties
- Reaction path calculations
- Heats of formation
- Dipole moments
- Ionization potentials
- UV/Vis spectra (CI)

## Key Strengths

### Speed:
- Semi-empirical efficiency
- Thousands of atoms feasible
- Fast screening
- Real-time geometry optimization
- Rapid conformational analysis

### AM1/PM3 Methods:
- Well-parametrized for organics
- Reasonable geometries
- Heat of formation accuracy
- Broad element coverage
- Proven performance

### Organic Chemistry Focus:
- Organic reactions
- Drug-like molecules
- Conformational analysis
- Biomolecular studies

### Practical Features:
- Transition states
- IRC following
- Thermochemistry
- Spectroscopy predictions

## Inputs & Outputs
- **Input formats**:
  - AMPAC input files
  - Molecular coordinates
  - Keyword-driven input
  
- **Output data types**:
  - Heats of formation
  - Optimized geometries
  - Frequencies
  - Properties

## Interfaces & Ecosystem
- **Standalone**: Complete package
- **Visualization**: Various molecular viewers
- **Successor**: MOPAC continuation

## Advanced Features

### AM1 Method:
- Core-core repulsion refinement
- Improved hydrogen bonding
- Better activation energies
- Broad parametrization

### PM3 Method:
- Reoptimized parameters
- Different functional form aspects
- Extended element coverage
- Complementary to AM1

### Reaction Calculations:
- Transition state optimization
- Intrinsic reaction coordinate
- Reaction paths
- Activation energies

### Spectroscopy:
- Vibrational frequencies
- IR intensities
- UV/Vis via CI
- Thermochemistry

## Performance Characteristics
- **Speed**: Orders of magnitude faster than ab initio
- **Accuracy**: Semi-empirical level (~5- 10 kcal/mol)
- **System size**: Thousands of atoms
- **Memory**: Minimal requirements

## Computational Cost
- **Energy**: Milliseconds
- **Optimization**: Seconds to minutes
- **Frequencies**: Minutes
- **Large systems**: Still fast
- **Typical**: Quick screening tool

## Limitations & Known Constraints
- **Accuracy**: Semi-empirical limitations
- **Elements**: Parametrization coverage
- **Exotic systems**: May fail outside training set
- **Superseded**: MOPAC is active successor
- **Availability**: Limited/historic

## Comparison with Other Codes
- **vs MOPAC**: MOPAC is active successor
- **vs xTB**: GFN-xTB more modern semi-empirical
- **vs PM7**: PM7 in MOPAC is evolved method
- **vs Ab initio**: Much faster, less accurate
- **Legacy**: Methods still widely used

## Application Areas

### Drug Discovery:
- Rapid screening
- Conformational analysis
- QSAR descriptors
- Lead optimization

### Organic Reactions:
- Mechanism exploration
- Transition state estimation
- Reaction energetics
- Selectivity prediction

### Large Molecules:
- Biomolecules
- Polymers
- Materials screening
- Supramolecular systems

### Education:
- Teaching quantum chemistry
- Demonstrating concepts
- Student projects
- Quick calculations

## Historical Context

### Development:
- 1980s: AM1 development (Dewar)
- 1989: PM3 introduction (Stewart)
- 1990s: Widespread adoption
- 2000s+: Continued in MOPAC

### Key Publications:
- Dewar AM1 paper (JACS 1985)
- Stewart PM3 paper (JCC 1989)
- Thousands of applications

## Community and Support
- Historic commercial product
- MOPAC as open continuation
- Extensive literature
- Methods standard in field
- Stewart continued development

## Verification & Sources
**Primary sources**:
1. Dewar et al., JACS 107, 3902 (1985) - AM1
2. Stewart, J. Comput. Chem. 10, 209 (1989) - PM3
3. MOPAC continuation: https://openmopac.net/
4. Extensive application literature

**Confidence**: VERIFIED (Historic)
- Status: Historic, succeeded by MOPAC
- Significance: AM1/PM3 methods foundational
- Impact: Thousands of applications
- Methods: Still widely used
