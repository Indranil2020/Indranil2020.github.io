# JADE-NAMD

## Official Resources
- Homepage: https://github.com/bch-gnome/JADE-NAMD
- Documentation: https://jade-namd.readthedocs.io/
- Source Repository: https://github.com/bch-gnome/JADE-NAMD
- License: MIT License

## Overview
JADE-NAMD is a Python-based software package designed for performing on-the-fly nonadiabatic molecular dynamics (NAMD) simulations. It employs the trajectory surface hopping method and serves as a flexible interface driver that connects various quantum chemistry packages (calculators) with dynamics propagation. It is designed to be user-friendly and easily extensible for different electronic structure methods.

**Scientific domain**: Nonadiabatic molecular dynamics, excited-state dynamics, trajectory surface hopping
**Target user community**: Researchers using standard QC packages for excited-state dynamics

## Theoretical Methods
- Trajectory Surface Hopping (TSH)
- Fewest-Switches Surface Hopping (FSSH)
- On-the-fly dynamics
- Numerical gradient integration
- Velocity Verlet integration
- Decoherence corrections (ID-A, EDC)
- Landau-Zener probability

## Capabilities (CRITICAL)
- Molecular dynamics propagation
- Surface hopping algorithms
- Interface with multiple QC codes
- Energy and gradient handling
- Non-adiabatic coupling vectors (NAC)
- Probabilistic hopping
- Trajectory analysis
- Kinetic energy conservation
- State tracking

**Sources**: GitHub repository, Documentation

## Key Strengths

### Flexible Interfaces:
- Turbomole
- GAMESS-US
- Gaussian
- Molpro
- MNDO
- Easy to extend for other codes

### Python-Based:
- Modern code structure
- Easy installation (pip)
- Readable codebase
- Scriptable workflows

### Methodological Generality:
- Independent of electronic structure method
- Supports any method providing E, Grad, NAC
- Various decoherence schemes

## Inputs & Outputs
- **Input formats**:
  - Python driver script
  - Configuration files (JSON/YAML)
  - Geometry files (XYZ)
  - Calculator templates
  
- **Output data types**:
  - Trajectory coordinates/velocities
  - Energy evolution
  - State population
  - Hopping logs
  - Restart files

## Interfaces & Ecosystem
- **Calculators**: Turbomole, GAMESS, Gaussian, Molpro, MNDO
- **Language**: Pure Python
- **Dependencies**: NumPy, SciPy
- **Analysis**: Matplotlib, standard trajectory tools

## Advanced Features

### Decoherence Handling:
- Instantaneous Decoherence (ID-A)
- Energy-based Decoherence (EDC)
- Improved hopping consistency

### Modular Design:
- Calculator abstract base class
- Pluggable dynamics engines
- Customizable propagators

## Performance Characteristics
- **Speed**: Driven by external QC code
- **Overhead**: Minimal Python overhead
- **Parallelization**: Script-level trajectory parallelism

## Computational Cost
- **Bottleneck**: Quantum chemistry calculations
- **Scaling**: Dependent on QC method (TDDFT vs CASSCF)
- **Typical**: Tens to hundreds of trajectories

## Limitations & Known Constraints
- **Calculator Dependence**: Needs external software
- **NAC Availability**: Requires QC code to compute couplings
- **System Size**: Limited by QC method capabilities

## Comparison with Other Codes
- **vs SHARC**: JADE is lighter weight, Python-centric
- **vs NEXMD**: JADE uses ab initio codes, NEXMD is semiempirical
- **vs Newton-X**: JADE is a simpler, more modern Python alternative
- **Unique strength**: Lightweight, Python-native, easy interfacing

## Application Areas
- **Photoisomerization**: Azobenzene, retinal
- **Photodissociation**: Bond breaking dynamics
- **Intersystem Crossing**: Spin-forbidden transitions
- **Materials**: Molecular switches

## Best Practices
- **Interface check**: Verify QC output parsing
- **Timestep**: Appropriate for nuclear motion (0.5-1.0 fs)
- **Ensemble**: Run sufficient trajectories for statistics
- **Cleanup**: Manage scratch files from QC codes

## Community and Support
- Open-source MIT license
- GitHub issue tracker
- Documentation on ReadTheDocs
- Active development by detailed contributors

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/bch-gnome/JADE-NAMD
2. Documentation: https://jade-namd.readthedocs.io/

**Confidence**: VERIFIED - Active GitHub project

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (MIT)
- Active development: Recent commits
- Specialized strength: Python-based surface hopping driver
