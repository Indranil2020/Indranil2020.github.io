# TDEP (Temperature Dependent Effective Potential)

## Official Resources
- Homepage: https://ollehellman.github.io/
- Documentation: https://ollehellman.github.io/page/documentation.html
- Source Repository: https://github.com/ollehellman/TDEP
- License: MIT License

## Overview
TDEP (Temperature Dependent Effective Potential) is software for extracting temperature-dependent force constants and studying anharmonic lattice dynamics using ab-initio molecular dynamics. Developed by Olle Hellman (Linköping University), TDEP uses temperature-dependent effective harmonic theory to capture anharmonic effects, making it powerful for materials with strong temperature dependence, soft phonon modes, and systems where perturbation theory fails.

**Scientific domain**: Temperature-dependent lattice dynamics, anharmonic phonons, thermal expansion  
**Target user community**: Phonon researchers, materials scientists studying temperature effects

## Theoretical Methods
- Temperature-dependent effective potential
- Effective harmonic theory  
- Force constant extraction from MD
- Self-consistent phonon theory
- Thermal expansion and free energy
- Phonon linewidths and lifetimes
- Renormalized phonons
- Grüneisen parameters

## Capabilities (CRITICAL)
- Temperature-dependent phonon spectra from AIMD
- Effective force constants extraction
- Thermal expansion coefficients
- Free energy vs temperature
- Phonon lifetimes and linewidths
- Soft mode stabilization
- Phase transition characterization
- Anharmonic renormalization
- VASP, QE, other DFT compatibility
- Handles strong anharmonicity

**Sources**: TDEP documentation, Phys. Rev. B 84, 180301(R) (2011); Phys. Rev. B 88, 144301 (2013)

## Key Strengths
- **Temperature-dependent**: True temperature effects from MD
- **AIMD-based**: Captures full anharmonicity
- **Soft modes**: Handles dynamical instabilities
- **Phase transitions**: Effective for structural transitions

## Inputs & Outputs
- **Input formats**: AIMD trajectories (VASP, QE), crystal structures, forces/positions
- **Output data types**: Temperature-dependent phonons, force constants, free energy, thermal expansion

## Interfaces & Ecosystem
- **VASP**: Primary interface
- **Quantum ESPRESSO**: Compatible
- **Fortran**: Core implementation
- **Python**: Post-processing tools
- **phonopy**: Visualization integration


## Advanced Features
- **Temperature-dependent effective potential**: True finite-temperature phonons
- **AIMD-based extraction**: Captures full anharmonicity from MD
- **Soft mode stabilization**: Handles imaginary phonon modes
- **Free energy calculations**: Thermodynamic properties vs temperature
- **Thermal expansion**: Grüneisen parameters and expansion coefficients
- **Phase transition detection**: Identifies structural instabilities
- **Phonon spectral functions**: Beyond harmonic approximation

## Performance Characteristics
- AIMD: Days (computationally expensive)
- TDEP processing: Minutes (fast)
- Overall: MD cost dominates

## Computational Cost
- DFT-MD: Dominant cost (days to weeks)
- TDEP extraction: Fast (minutes to hours)
- Separate MD run needed per temperature

## Limitations & Known Constraints
- **Requires expensive AIMD**: Many MD steps needed
- **MD convergence critical**: Sufficient sampling required
- **System size limitations**: MD supercell constraints
- **Learning curve**: Moderate
- **Temperature scanning**: Each T requires separate MD

## Comparison with Other Codes
- **vs SSCHA**: Both handle strong anharmonicity; different methodologies
- **vs perturbative phonons**: TDEP for strong temperature dependence
- **Unique strength**: Temperature-dependent effective potential from MD

## Application Areas
- Temperature-dependent phonon spectroscopy
- Soft phonon mode materials
- Structural phase transitions
- Thermal expansion studies
- Thermoelectric materials
- High-temperature phonon physics

## Best Practices
- Sufficient MD statistics (1000+ steps minimum)
- Converge supercell size
- Multiple temperatures for phase diagrams
- Validate against experimental phonon data
- Check force constant convergence

## Community and Support
- Open-source (MIT license)
- GitHub repository
- Documentation website
- Active development
- Growing user base

## Educational Resources
- Comprehensive documentation
- Tutorial examples
- Publications with methodology
- Example calculations

## Development
- Olle Hellman (Linköping University, Sweden)
- Active development
- Regular updates
- Well-maintained

## Research Impact
TDEP enables accurate temperature-dependent phonon calculations from AIMD, crucial for materials with strong anharmonicity, soft modes, and temperature-driven phase transitions.

## Verification & Sources
**Primary sources**:
1. Homepage: https://ollehellman.github.io/
2. Documentation: https://ollehellman.github.io/page/documentation.html
3. GitHub: https://github.com/ollehellman/TDEP
4. Publications: Phys. Rev. B 84, 180301(R) (2011); Phys. Rev. B 88, 144301 (2013)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE and ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub, MIT license)
- Development: ACTIVE (Linköping University)
- Publications: PEER-REVIEWED
- Applications: Temperature-dependent phonons, MD-based effective potential, strong anharmonicity, soft modes, phase transitions, thermal expansion, production quality
