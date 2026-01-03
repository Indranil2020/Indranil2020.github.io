# ALM (Anharmonic Lattice Model - Force Constant Extraction)

## Official Resources
- Homepage: https://github.com/ttadano/ALM (part of ALAMODE project)
- Documentation: https://alamode.readthedocs.io/en/latest/alm.html
- Source Repository: https://github.com/ttadano/alamode
- License: MIT License

## Overview
ALM is the force constant extraction module of the ALAMODE software suite. It extracts harmonic and anharmonic interatomic force constants from first-principles displacement-force datasets using advanced fitting techniques including compressive sensing. ALM can extract 2nd, 3rd, 4th, and higher-order force constants efficiently, making it a key tool for setting up anharmonic lattice dynamics calculations.

**Scientific domain**: Force constant extraction, anharmonic lattice dynamics  
**Target user community**: Researchers studying phonons and thermal transport

## Theoretical Methods
- Least-squares fitting
- Compressive sensing (LASSO, elastic net)
- Symmetry constraints
- Cutoff-based cluster expansion
- Harmonic and anharmonic IFCs (2nd through 6th order)
- Regularization techniques
- Cross-validation

## Capabilities (CRITICAL)
- Extract 2nd order (harmonic) force constants
- Extract 3rd, 4th, 5th, 6th order anharmonic force constants
- Compressive sensing for data efficiency
- Automatic symmetry implementation
- Compatible with any DFT code (via force-displacement data)
- Efficient algorithms for large systems
- Python and C++ interfaces
- Integration with ALAMODE for thermal conductivity

**Sources**: ALAMODE/ALM documentation, J. Phys.: Condens. Matter 26, 225402 (2014)

## Key Strengths
- **High-order IFCs**: Up to 6th order anharmonicity
- **Compressive sensing**: Reduces required calculations
- **Flexibility**: Works with any DFT code
- **Symmetry**: Automatic symmetry constraints

## Inputs & Outputs
- **Input formats**:
  - Displacement patterns
  - Forces from DFT calculations
  - Crystal structure
  - Symmetry information
  
- **Output data types**:
  - Harmonic force constants (phonopy format)
  - Anharmonic force constants (various orders)
  - ALAMODE format for thermal conductivity

## Interfaces & Ecosystem
- **ALAMODE**: Core component for full thermal transport workflow
- **phonopy**: Export harmonic force constants
- **DFT codes**: Any code (via displacement-force data)
- **Python/C++**: Both interfaces available

## Workflow and Usage

### Extract Harmonic Force Constants:
```bash
# Generate displacements
alm -p input.in --suggest

# After DFT calculations
alm input.in
```

### Extract Anharmonic (3rd order):
```
&general
  PREFIX = silicon
  MODE = optimize
  NAT = 8; NKD = 1
  KD = Si
/

&interaction
  NORDER = 2  # 2=3rd order
/

&optimize
  LMODEL = enet  # elastic net regression
/
```

## Advanced Features
- Multiple regression methods (least-squares, LASSO, elastic net)
- Cross-validation for model selection
- Sparse optimization techniques
- High-order force constants for extreme anharmonicity

## Performance Characteristics
- Fast fitting algorithms
- Handles large displacement datasets
- Efficient memory usage
- Parallelization support

## Computational Cost
- DFT calculations: Dominant cost
- ALM fitting: Fast (seconds to minutes)
- Scales well with system size

## Limitations & Known Constraints
- **Requires DFT input**: Not a DFT code itself
- **Displacement generation**: Manual or external tools
- **Learning curve**: Moderate
- **High-order IFCs**: Require many calculations

## Comparison with Other Codes
- **Part of ALAMODE**: Integrated force constant extraction module
- **vs hiPhive**: Similar compressive sensing; ALM more established
- **vs phono3py direct fit**: ALM more flexible algorithms

## Application Areas
- Anharmonic phonon calculations
- Thermal conductivity (via ALAMODE)
- Force constant database generation
- High-throughput phonon studies

## Best Practices
- Use compressive sensing to reduce calculations
- Systematic convergence of cutoffs
- Cross-validation to prevent overfitting
- Test different regression methods

## Community and Support
- Part of ALAMODE project (MIT license)
- GitHub repository
- ALAMODE documentation
- Active development
- User support via issues

## Development
- Terumasa Tadano (Tohoku University)
- Part of ALAMODE development
- Regular updates
- Well-maintained

## Research Impact
ALM provides efficient extraction of anharmonic force constants using advanced regression techniques, enabling accurate thermal transport calculations with reduced computational cost.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ttadano/alamode
2. Documentation: https://alamode.readthedocs.io/en/latest/alm.html
3. Publication: J. Phys.: Condens. Matter 26, 225402 (2014)

**Confidence**: VERIFIED - Part of ALAMODE

**Verification status**: ✅ VERIFIED
- Part of ALAMODE suite (MIT license)
- Documentation: COMPREHENSIVE
- Development: ACTIVE (Tohoku)
- Applications: Force constant extraction, compressive sensing, harmonic and anharmonic IFCs, integration with ALAMODE thermal transport
