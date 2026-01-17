# SHARC (Surface Hopping including ARbitrary Couplings)

## Official Resources
- Homepage: https://sharc-md.org/
- Documentation: https://sharc-md.org/?page_id=50
- Source Repository: https://github.com/sharc-md/sharc
- License: GNU General Public License v3.0

## Overview
SHARC is a comprehensive ab initio molecular dynamics software suite for excited-state dynamics simulations using trajectory surface hopping. It enables the study of photochemical and photophysical processes including internal conversion, intersystem crossing, and photodissociation. SHARC interfaces with major quantum chemistry codes to obtain electronic structure data and supports various types of couplings including spin-orbit and non-adiabatic couplings.

**Scientific domain**: Photochemistry, excited-state dynamics, nonadiabatic processes, intersystem crossing
**Target user community**: Researchers studying light-induced chemical reactions, photophysics, and ultrafast dynamics

## Theoretical Methods
- Trajectory Surface Hopping (TSH)
- Fewest-Switches Surface Hopping (FSSH)
- Landau-Zener surface hopping
- Non-adiabatic coupling vectors
- Spin-orbit couplings (SOC)
- Decoherence corrections
- Velocity rescaling and adjustment
- SHARC propagator (diagonal representation)
- MCH propagator (molecular Coulomb Hamiltonian)

## Capabilities (CRITICAL)
- Excited-state molecular dynamics
- Non-adiabatic dynamics simulations
- Intersystem crossing (ISC) dynamics
- Internal conversion pathways
- Spin-orbit coupled dynamics
- Wigner sampling for initial conditions
- Newton-X interface for initial conditions
- Trajectory analysis and plotting
- Ensemble averaging and statistics
- Time-resolved property propagation
- Dipole moments and transition properties
- Natural transition orbitals
- Population dynamics

**Sources**: Official SHARC documentation, GitHub repository, J. Chem. Theory Comput. 2019

## Key Strengths

### Arbitrary Couplings:
- Non-adiabatic couplings
- Spin-orbit couplings
- Dipole couplings (laser fields)
- Simultaneous treatment of all couplings
- Flexible coupling definitions

### Extensive QC Interfaces:
- OpenMolcas (CASSCF/CASPT2/MS-CASPT2)
- ORCA (TDDFT, MRCI)
- TURBOMOLE (TDDFT, CC2, ADC(2))
- ADF (TDDFT)
- Gaussian (TDDFT)
- Columbus (MRCI)
- BAGEL (CASPT2)

### Analysis Tools:
- Trajectory visualization
- Population analysis
- Reaction pathway analysis
- Statistical convergence
- Property time-evolution

### Machine Learning Integration:
- SchNarc interface
- ML potential energy surfaces
- Neural network dynamics
- Accelerated sampling

## Inputs & Outputs
- **Input formats**:
  - SHARC input files
  - Geometry files (XYZ)
  - Wigner ensemble distributions
  - QC interface templates
  
- **Output data types**:
  - Trajectory files
  - Population dynamics
  - Geometries and velocities
  - Electronic state data
  - Hopping statistics
  - Property evolution

## Interfaces & Ecosystem
- **QC programs**: OpenMolcas, ORCA, TURBOMOLE, ADF, Gaussian, Columbus, BAGEL
- **Visualization**: VMD, PyMOL, trajectory viewers
- **ML integration**: SchNarc (neural network potentials)
- **Initial conditions**: Newton-X interface, Wigner sampling
- **Post-processing**: Python analysis scripts, matplotlib

## Advanced Features

### SHARC Propagator:
- Diagonal representation of electronic states
- Proper treatment of couplings
- Energy conservation
- Momentum adjustment schemes

### Laser Field Dynamics:
- Time-dependent electric fields
- Dipole coupling terms
- Photoexcitation modeling
- Pulse shape control

### Adaptive Sampling:
- Efficient trajectory selection
- Importance sampling
- Enhanced statistics

## Performance Characteristics
- **Speed**: Determined by QC interface
- **Accuracy**: Based on underlying electronic structure method
- **System size**: Limited by QC method (typically <500 atoms)
- **Parallelization**: Embarrassingly parallel trajectory ensemble

## Computational Cost
- **Dynamics**: Overhead minimal compared to QC
- **Bottleneck**: Electronic structure calculations
- **Typical**: Hundreds to thousands of trajectories
- **Resources**: HPC cluster recommended

## Limitations & Known Constraints
- **Electronic structure**: Requires external QC program
- **System size**: Limited by QC method
- **Classical nuclei**: Ignores nuclear quantum effects
- **Decoherence**: Model-dependent corrections
- **Basis set**: Must be consistent across trajectory

## Comparison with Other Codes
- **vs Newton-X**: SHARC handles arbitrary couplings, both TSH
- **vs NEXMD**: SHARC more general interfaces, NEXMD semiempirical focus
- **vs CPMD/CP2K**: SHARC specialized for hopping, general codes have Ehrenfest
- **Unique strength**: Spin-orbit couplings, arbitrary coupling treatment, extensive interfaces

## Application Areas

### Photochemistry:
- Bond photodissociation
- Ring opening/closing reactions
- Proton/hydrogen transfer
- Isomerization dynamics

### Photophysics:
- Intersystem crossing rates
- Phosphorescence mechanisms
- Triplet state dynamics
- Heavy-atom effects

### Biological Systems:
- DNA photodamage
- Retinal photoisomerization
- Photoreceptor proteins
- Chromophore dynamics

## Best Practices

### Method Selection:
- CASPT2/MS-CASPT2 for accuracy
- TDDFT for larger systems
- Careful active space selection
- Validate with static calculations

### Ensemble Size:
- Statistical convergence analysis
- Typically 100-500 trajectories
- Check property convergence
- Consider trajectory weights

### Analysis:
- Population dynamics
- Branching ratios
- Time constants
- Mechanism identification

## Community and Support
- Open-source GPL v3
- Active development (González and Truhlar groups)
- Tutorial materials available
- Regular workshops
- Published methodology papers

## Verification & Sources
**Primary sources**:
1. Official website: https://sharc-md.org/
2. GitHub repository: https://github.com/sharc-md/sharc
3. S. Mai et al., WIREs Comput. Mol. Sci. 8, e1370 (2018)
4. M. Richter et al., J. Chem. Theory Comput. 7, 1253 (2011)

**Secondary sources**:
1. SHARC tutorials and documentation
2. Published applications
3. Workshop materials

**Confidence**: VERIFIED - Active open-source project

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub, GPL v3)
- Community support: Active (developers, tutorials)
- Academic citations: >500
- Active development: Regular releases
- Specialized strength: Arbitrary couplings, extensive QC interfaces, spin-orbit dynamics
