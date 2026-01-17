# PyFLOSIC

## Official Resources
- Homepage: https://github.com/pyflosic/pyflosic
- Documentation: https://pyflosic.readthedocs.io/
- Source Repository: https://github.com/pyflosic/pyflosic
- PyPI: https://pypi.org/project/pyflosic/
- License: Apache License 2.0

## Overview
PyFLOSIC is a Python implementation of the Fermi-Löwdin Orbital Self-Interaction Correction (FLO-SIC) method built on top of PySCF. It provides an accessible interface for performing self-interaction corrected DFT calculations, improving orbital energies and other properties affected by self-interaction error.

**Scientific domain**: Molecules, self-interaction correction, accurate orbital energies  
**Target user community**: PySCF users needing self-interaction correction, researchers studying orbital energetics

## Theoretical Methods
- Perdew-Zunger Self-Interaction Correction (PZ-SIC)
- Fermi-Löwdin Orbital transformation
- Fermi Orbital Descriptors (FODs)
- Integration with PySCF functionals
- LDA, GGA, meta-GGA support
- SCF with SIC

## Capabilities (CRITICAL)
- Self-interaction corrected DFT
- Improved orbital energies
- FOD optimization
- PySCF integration
- Ionization potentials
- Electron affinities
- Multiple functional support
- Visualization of FODs
- Restart capabilities

**Sources**: GitHub repository, ReadTheDocs

## Key Strengths

### PySCF Foundation:
- Built on powerful PySCF
- Inherit PySCF capabilities
- Python ecosystem integration
- Extensive basis sets

### FOD Tools:
- FOD initialization
- Automatic placement
- Optimization algorithms
- Visualization

### Accessibility:
- pip installable
- Python interface
- Well documented
- Examples provided

## Inputs & Outputs
- **Input formats**:
  - PySCF Mole objects
  - FOD files
  - Python API
  
- **Output data types**:
  - SIC energies
  - Orbital energies
  - Optimized FODs
  - Properties

## Interfaces & Ecosystem
- **PySCF integration**:
  - Direct Mole/SCF use
  - Basis set support
  - Functional library
  
- **Visualization**:
  - FOD plotting
  - Orbital visualization
  - Integration with viewers

## Advanced Features

### FOD Optimization:
- Gradient-based
- Automatic initialization
- Constrained optimization
- Multiple algorithms

### Restart:
- Save/load FODs
- Checkpoint calculations
- Continue optimization

## Performance Characteristics
- **Speed**: PySCF backend efficiency
- **Accuracy**: Improved over standard DFT
- **System size**: Small to medium molecules
- **Memory**: PySCF requirements

## Limitations & Known Constraints
- **System size**: Molecular focus
- **Periodicity**: Not supported
- **Scaling**: Additional SIC overhead
- **FOD quality**: Initialization-dependent

## Comparison with Other Codes
- **vs FLOSIC**: PyFLOSIC Python/PySCF, FLOSIC Fortran
- **vs Standard DFT**: Self-interaction corrected
- **Unique strength**: Python accessibility, PySCF integration

## Application Areas

### Orbital Energetics:
- HOMO-LUMO gaps
- Ionization potentials
- Spectroscopy predictions
- Photoelectron spectra
- Koopmans' theorem validation

### Charge Transfer:
- Localized states
- Electron transfer
- Redox chemistry
- Donor-acceptor systems
- Charge-transfer excitations

### Barrier Heights:
- Reaction barriers
- SIE correction improves barriers
- Kinetics predictions
- Transition state energies

### Molecular Properties:
- Dipole moments
- Polarizabilities
- Magnetic properties
- Electronegativity scales

## Best Practices

### FOD Initialization:
- Start with core-like positions
- Use molecular symmetry
- Chemical intuition for lone pairs
- Test multiple starting guesses

### Convergence:
- Monitor SIC energy decrease
- Check FOD position stability
- Use tight convergence criteria
- Verify against analytical gradients

### Functional Selection:
- LDA-SIC well characterized
- GGA-SIC (PBE-SIC) commonly used
- Meta-GGA-SIC available
- Document functional choice

### Validation:
- Compare with experimental IPs
- Check against higher-level theory
- Test on small benchmark systems
- Monitor total energy consistency

## Computational Cost
- **SIC overhead**: 2-4x over standard DFT per iteration
- **FOD optimization**: Adds 5-20 additional optimization cycles
- **Per-orbital**: SIC applied to each occupied orbital
- **Scaling**: Cubic O(N³) with number of electrons
- **Memory**: PySCF-level requirements
- **Typical runs**: Minutes for small molecules, hours for medium

## Community and Support
- Open source Apache 2.0
- GitHub development
- ReadTheDocs
- PyPI distribution
- Active maintenance

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/pyflosic/pyflosic
2. ReadTheDocs: https://pyflosic.readthedocs.io/
3. PyPI: https://pypi.org/project/pyflosic/

**Confidence**: VERIFIED - Active Python package

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, Apache 2.0)
- Package: PyPI
- Documentation: ReadTheDocs
- Specialty: Self-interaction correction via PySCF
