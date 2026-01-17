# NEXMD (Nonadiabatic EXcited-state Molecular Dynamics)

## Official Resources
- Homepage: https://github.com/lanl/NEXMD
- Documentation: https://nexmd.github.io/
- Source Repository: https://github.com/lanl/NEXMD
- License: BSD 3-Clause License

## Overview
NEXMD is a software package developed at Los Alamos National Laboratory for simulating photoinduced adiabatic and non-adiabatic excited-state molecular dynamics. It uses semiempirical quantum chemistry methods (CEO package) with Tully's fewest-switches surface hopping algorithm, making it efficient for studying large conjugated systems like chromophores and polymers. Written in Fortran 90 with Python scripts for parallel execution.

**Scientific domain**: Photochemistry in large molecules, conjugated systems, chromophore dynamics, exciton dynamics
**Target user community**: Researchers studying organic chromophores, conjugated polymers, and large molecular photophysics

## Theoretical Methods
- Semiempirical methods (AM1, PM3, PM6, PM7)
- Collective Electronic Oscillator (CEO) approach
- Configuration Interaction Singles (CIS)
- Tully's Fewest-Switches Surface Hopping (FSSH)
- Non-adiabatic coupling calculations
- Trivial crossing detection
- Decoherence corrections (EDC, AFSSH)
- Velocity rescaling schemes

## Capabilities (CRITICAL)
- Ground and excited-state dynamics
- Non-adiabatic transitions between states
- Large molecular systems (100s of atoms)
- Exciton dynamics in conjugated systems
- Hot carrier relaxation
- Energy transfer simulations
- Absorption spectra simulation
- Time-resolved properties
- Parallel trajectory execution
- Ensemble averaging

**Sources**: Official GitHub documentation, LANL publications

## Key Strengths

### Semiempirical Efficiency:
- Fast electronic structure
- Large systems (>100 atoms)
- Many excited states feasible
- Long timescale dynamics

### CEO Methodology:
- Collective modes description
- Efficient gradient computation
- Multi-state treatment
- Excited-state forces

### LANL Development:
- National lab support
- Active maintenance
- Regular updates
- Scientific validation

### Conjugated Systems:
- Polymer dynamics
- Organic chromophores
- Exciton migration
- Energy transfer

## Inputs & Outputs
- **Input formats**:
  - NEXMD input files
  - Geometry files (XYZ)
  - Parameter files
  - Python driver scripts
  
- **Output data types**:
  - Trajectory files
  - Population dynamics
  - Energy files
  - Coupling data
  - Statistical output

## Interfaces & Ecosystem
- **Electronic structure**: Internal CEO (semiempirical)
- **Scripting**: Python for job management
- **Parallelization**: Multiple trajectory parallelism
- **Analysis**: Built-in analysis tools

## Advanced Features

### Trivial Crossing Detection:
- Automatic state relabeling
- Diabatic following
- Coupling analysis
- State tracking

### Decoherence Methods:
- Energy-based decoherence (EDC)
- Augmented FSSH (AFSSH)
- Wavefunction collapse schemes
- Physical decoherence treatment

### Large-Scale Dynamics:
- Efficient for 100+ atoms
- Many trajectories feasible
- Long timescales (ps)
- Statistical averaging

## Performance Characteristics
- **Speed**: Fast (semiempirical)
- **Accuracy**: Good for trends
- **System size**: 100s of atoms
- **Parallelization**: Trajectory-level

## Computational Cost
- **Per step**: Milliseconds (semiempirical)
- **Full trajectory**: Minutes to hours
- **Ensemble**: Highly parallel
- **Typical**: 100-500 trajectories

## Limitations & Known Constraints
- **Accuracy**: Semiempirical limitations
- **Parametrization**: Requires validated parameters
- **Heavy atoms**: Limited treatment
- **Spin-orbit**: Not included
- **Quantum nuclei**: Classical only

## Comparison with Other Codes
- **vs SHARC/Newton-X**: NEXMD faster but less accurate
- **vs DFTBaby**: Similar semiempirical approach
- **vs Ab initio codes**: NEXMD faster, less accurate
- **Unique strength**: Large conjugated systems, efficiency, LANL support

## Application Areas

### Conjugated Polymers:
- Polythiophenes
- PPV derivatives
- Donor-acceptor polymers
- Exciton dynamics

### Organic Chromophores:
- Dye molecules
- Photosensitizers
- OLED materials
- Photovoltaic materials

### Energy Transfer:
- FRET dynamics
- Exciton migration
- Hot carrier relaxation
- Charge separation

### Biological Chromophores:
- Chlorophylls
- Carotenoids
- Flavins
- Photoactive proteins

## Best Practices

### Parameter Validation:
- Benchmark against ab initio
- Check excited-state ordering
- Validate geometries
- Compare spectroscopy

### Trajectory Management:
- Sufficient ensemble size
- Convergence checking
- Statistical analysis
- Error estimation

### System Setup:
- Proper initial conditions
- Adequate equilibration
- Careful state selection
- Documentation

## Community and Support
- Open-source BSD license
- LANL development team
- GitHub repository
- Documentation available
- Active maintenance

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/lanl/NEXMD
2. Documentation: https://nexmd.github.io/
3. S. Tretiak et al., publications on CEO methodology
4. A. F. Fidler et al., J. Phys. Chem. Lett. 2013

**Secondary sources**:
1. LANL software registry
2. Published applications
3. Conference presentations

**Confidence**: VERIFIED - LANL open-source

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE
- Source code: OPEN (BSD 3-Clause)
- Community support: LANL maintained
- Active development: Yes
- Specialized strength**: Large conjugated systems, semiempirical efficiency, exciton dynamics
