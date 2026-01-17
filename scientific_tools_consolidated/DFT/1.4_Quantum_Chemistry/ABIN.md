# ABIN

## Official Resources
- Homepage: https://github.com/PHOTOX/ABIN
- Documentation: https://github.com/PHOTOX/ABIN/wiki
- Source Repository: https://github.com/PHOTOX/ABIN
- License: GNU General Public License v3.0

## Overview
ABIN (Ab Initio Born-oppenheimer Nuclear dynamics) is a multipurpose ab initio molecular dynamics program. It is designed to perform ab initio MD and model nuclear quantum effects, interfacing with external electronic structure programs like ORCA and TeraChem for forces and energies.

**Scientific domain**: Ab initio molecular dynamics, nuclear quantum effects  
**Target user community**: Researchers studying molecular dynamics with quantum nuclear effects

## Theoretical Methods
- Ab initio molecular dynamics (AIMD)
- Path Integral Molecular Dynamics (PIMD)
- Ring Polymer Molecular Dynamics (RPMD)
- Surface hopping dynamics
- Centroid molecular dynamics
- Multiple electronic structure backends

## Capabilities (CRITICAL)
- Born-Oppenheimer MD
- Path integral nuclear quantization
- Surface hopping for nonadiabatic dynamics
- Multiple replica propagation
- Ring polymer methods
- Interface to ORCA, TeraChem, Gaussian
- NVE, NVT, NPT ensembles
- Thermostat implementations
- Trajectory analysis

## Key Strengths

### Nuclear Quantum Effects:
- Path integrals
- Ring polymer MD
- Centroid dynamics
- Quantum tunneling
- Zero-point energy

### Interface Architecture:
- Shell script interface
- Multiple QC backends
- Easy code swapping
- Flexible input

### Nonadiabatic Dynamics:
- Surface hopping
- Multiple states
- Excited state dynamics
- Photochemistry

### MD Capabilities:
- Standard integrators
- Thermostats
- Barostats
- Trajectory output

## Inputs & Outputs
- **Input formats**:
  - ABIN input files
  - Coordinates (XYZ)
  - Velocities
  
- **Output data types**:
  - Trajectories
  - Energies/forces
  - PIMD observables
  - Statistical properties

## Interfaces & Ecosystem
- **QC backends**: ORCA, TeraChem, Gaussian, Molpro
- **Analysis**: Trajectory tools
- **Visualization**: Standard MD formats

## Advanced Features

### Path Integrals:
- Bead propagation
- Staging coordinates
- PILE thermostat
- Convergence with beads

### Surface Hopping:
- Tully's FSSH
- Multiple states
- Decoherence corrections
- Hopping algorithms

### Replica Methods:
- Multiple trajectory
- Parallel execution
- Ensemble averaging
- Uncertainty quantification

## Performance Characteristics
- **Speed**: QC-limited
- **Accuracy**: Backend accuracy
- **System size**: Moderate (QC limited)
- **Parallelization**: Replica parallel

## Computational Cost
- **Classical AIMD**: QC cost per step
- **PIMD**: Beads × QC cost
- **Surface hopping**: Multiple states
- **Typical**: QC is bottleneck

## Limitations & Known Constraints
- **Electronic structure**: External dependency
- **Large systems**: QC limitations
- **Documentation**: Research-focused
- **Setup**: Interface configuration needed

## Comparison with Other Codes
- **vs i-PI**: Both PIMD; different interfaces
- **vs CP2K**: ABIN lighter, backend-agnostic
- **vs Newton-X**: Both nonadiabatic; different focus
- **Unique strength**: PIMD + QC interfaces

## Application Areas

### Photochemistry:
- Excited state dynamics
- Photodissociation
- Internal conversion
- Intersystem crossing

### Nuclear Quantum Effects:
- Hydrogen transfer
- Tunneling reactions
- Isotope effects
- Light atom dynamics

### Condensed Phase:
- Solutions
- Interfaces
- Quantum solvent effects

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/PHOTOX/ABIN
2. PHOTOX group (Charles University)
3. Related publications

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Documentation: Wiki
- Active development: Yes
