# FLOSIC (Fermi-Löwdin Orbital Self-Interaction Correction)

## Official Resources
- Homepage: https://www.flosic.org/
- Documentation: https://flosic.utep.edu/
- Source Repository: https://github.com/FLOSIC
- License: Open Source (DOE funded)

## Overview
FLOSIC is an electronic structure code that implements the Fermi-Löwdin Orbital Self-Interaction Correction (FLO-SIC) method to address self-interaction errors in standard DFT calculations. Built upon NRLMOL, it uses Gaussian orbitals and provides improved predictions for orbital energies, ionization potentials, and electron affinities.

**Scientific domain**: Molecules, reaction barriers, redox chemistry, orbital energetics  
**Target user community**: Researchers needing accurate orbital energies, ionization potentials, and self-interaction-free DFT

## Theoretical Methods
- Density Functional Theory (DFT)
- Fermi-Löwdin Orbital Self-Interaction Correction (FLO-SIC)
- Perdew-Zunger SIC formulation
- Gaussian orbital basis sets
- LDA and GGA exchange-correlation functionals
- Fermi orbital descriptors (FODs)
- Self-consistent SIC implementation

## Capabilities (CRITICAL)
- Ground-state electronic structure with SIC
- Self-interaction corrected total energies
- Improved orbital energies
- Accurate ionization potentials
- Electron affinities
- Reaction barriers
- Charge transfer states
- Redox potentials
- Fermi orbital descriptor optimization
- Massively parallel calculations

**Sources**: FLOSIC Center, UTEP, DOE Computational Chemistry Program

## Key Strengths

### Self-Interaction Correction:
- Removes one-electron self-interaction errors
- Improved orbital energetics
- Better HOMO-LUMO gaps
- More physical description of electrons

### Fermi-Löwdin Orbitals:
- Transformation from KS orbitals
- Orthogonalized Fermi orbitals
- Unique SIC energy evaluation
- Physical interpretation

### FOD Optimization:
- Electronic geometry optimization
- Fermi orbital descriptors in 3D
- Minimizes SIC energy
- Automatic or manual placement

### Broad Applicability:
- Molecules and clusters
- Charged species
- Transition states
- Radical chemistry

## Inputs & Outputs
- **Input formats**:
  - CLUSTER file (geometry)
  - FRMORB file (FOD positions)
  - Basis set specifications
  - Control parameters
  
- **Output data types**:
  - SIC-corrected energies
  - Orbital energies
  - Optimized FOD positions
  - Forces
  - Charge analysis

## Interfaces & Ecosystem
- **NRLMOL base**:
  - Built on Naval Research Lab code
  - UTEP modifications for FLO-SIC
  - Fortran implementation
  
- **PyFLOSIC**:
  - Python implementation
  - PySCF integration
  - Simplified interface

## Advanced Features

### FOD Optimization:
- Gradient-based optimization
- Automatic initial guesses
- Constrained optimization
- Multiple starting points

### Parallel Implementation:
- MPI parallelization
- Massively parallel scaling
- Large system capability
- HPC-ready

### Multiple SIC Schemes:
- Standard PZ-SIC
- Scaled SIC variants
- Functional-specific corrections

### Property Calculations:
- Ionization potentials from orbital energies
- Electron affinities
- Koopman's-like behavior
- Photoemission predictions

## Performance Characteristics
- **Speed**: Moderate (SIC adds overhead)
- **Accuracy**: Improved over standard DFT for many properties
- **System size**: Small to medium molecules
- **Memory**: Standard Gaussian code requirements
- **Parallelization**: Excellent MPI scaling

## Computational Cost
- **SIC overhead**: 2-4x standard DFT
- **FOD optimization**: Additional iterations
- **Scaling**: Cubic with system size
- **Typical**: Hours for medium molecules

## Limitations & Known Constraints
- **System size**: Best for molecules (< 100-200 atoms)
- **FOD initialization**: Requires good starting positions
- **Periodicity**: Not periodic (molecular focus)
- **Functionals**: Not all functionals tested
- **Complexity**: Additional SIC parameters

## Comparison with Other Codes
- **vs Standard DFT**: FLOSIC corrects self-interaction
- **vs Hybrid DFT**: Different approach to exchange error
- **vs GW**: Both improve orbital energies
- **Unique strength**: Systematic self-interaction correction, FOD formalism

## Application Areas

### Redox Chemistry:
- Electron transfer reactions
- Redox potentials
- Oxidation states
- Battery materials

### Radical Chemistry:
- Open-shell systems
- Radical stability
- Spin contamination reduction
- Transition metal centers

### Ionization Potentials:
- Photoemission prediction
- HOMO energies
- Core ionization
- Valence shell

### Reaction Barriers:
- Transition state energetics
- SIE affects barriers
- Better kinetics predictions

## Best Practices

### FOD Placement:
- Chemical intuition helps
- Start with atomic cores
- Optimize thoroughly
- Check for local minima

### Functional Choice:
- Test with parent functional first
- LDA-SIC well characterized
- GGA-SIC available
- Document functional used

### Convergence:
- Monitor SIC energy convergence
- Check FOD movement
- Verify final positions reasonable

## Community and Support
- FLOSIC Center (UTEP, Central Michigan, Others)
- DOE Computational Chemistry
- GitHub repositories
- Published methodology
- Active development

## Verification & Sources
**Primary sources**:
1. FLOSIC website: https://www.flosic.org/
2. GitHub: https://github.com/FLOSIC
3. UTEP documentation: https://flosic.utep.edu/
4. M. R. Pederson et al., J. Chem. Phys. publications

**Confidence**: VERIFIED - DOE-funded, active development

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Academic use: Growing community
- Documentation: Available
- Active development: Regular updates
- Specialty: Self-interaction correction, orbital energetics
