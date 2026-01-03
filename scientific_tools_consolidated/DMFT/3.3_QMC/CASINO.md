# CASINO (Cambridge Stochastic Investigation of Novel Observables)

## Official Resources
- Homepage: https://vallico.net/casino/
- Documentation: https://vallico.net/casino/documentation/
- Source Repository: Available to registered users
- License: Academic license (free for academics)

## Overview
CASINO is a mature and feature-rich quantum Monte Carlo (QMC) code for electronic structure calculations, developed at the University of Cambridge and maintained by an international collaboration. Known for its extensive capabilities and rigorous implementation of QMC methods, CASINO provides VMC, DMC, and related techniques with numerous advanced features for molecules, solids, and surfaces. The code is widely used in the QMC community and represents decades of methodological development.

**Scientific domain**: Quantum Monte Carlo, electronic structure, many-body physics  
**Target user community**: QMC specialists, computational chemists, materials scientists

## Theoretical Methods
- Variational Monte Carlo (VMC)
- Diffusion Monte Carlo (DMC)
- Reptation QMC
- Fixed-node and released-node DMC
- Slater-Jastrow wavefunctions
- Geminal wavefunctions
- Backflow transformations
- Pairing wavefunctions
- Multi-reference expansions

## Capabilities (CRITICAL)
**Category**: Academic QMC code
- VMC and DMC methods
- Molecules and solids
- Periodic and finite systems
- Advanced trial wavefunctions
- Geminal and pairing functions
- Backflow correlations
- Excited states
- Forces and geometry optimization
- Polarization calculations
- Wavefunction analysis
- Pseudopotentials and all-electron
- MPI parallelization
- Production quality

**Sources**: Official website, documentation, publications

## Key Strengths

### Feature-Rich:
- Extensive wavefunction types
- Advanced correlation functions
- Many observables
- Comprehensive tools
- Research capabilities

### Mature Code:
- Decades of development
- Well-tested
- Extensive validation
- Production quality
- Large user base

### Methodological Depth:
- Advanced QMC methods
- Released-node DMC
- Geminal functions
- Pairing wavefunctions
- Cutting-edge features

### Community Standard:
- Widely cited
- Benchmark calculations
- Method development
- Educational use
- Publications

## Inputs & Outputs
- **Input formats**:
  - CASINO input files
  - DFT trial wavefunctions (CASTEP, CRYSTAL, Gaussian, etc.)
  - gwfn.data (wavefunction)
  - Pseudopotentials
  
- **Output data types**:
  - Energies and forces
  - Structural properties
  - Excited states
  - Expectation values
  - Density matrices
  - Analysis data

## Interfaces & Ecosystem

### DFT Codes:
- CASTEP (native)
- CRYSTAL
- Gaussian
- GAMESS-US
- Turbomole
- Many others via converters

### Tools:
- runqmc (job manager)
- casinohelp
- Analysis utilities
- Plotting tools

## Workflow and Usage

### Installation:
```bash
# Register and download from website
# Extract and compile
tar -xzf CASINO*.tar.gz
cd CASINO*
make
```

### DFT Trial Wavefunction:
```bash
# CASTEP example
castep.serial material
# Produces .check file
# Convert to CASINO format
```

### VMC Calculation:
```bash
# Edit input file
cat > input << EOF
runtype : vmc_opt
neu : 4
ned : 4
atom_basis_type : gaussian
EOF

# Run CASINO
casino > vmc.out
```

### DMC Calculation:
```bash
# DMC input
cat > input << EOF
runtype : vmc_dmc
neu : 4
ned : 4
dmc_target_weight : 1000
dmc_equil_nstep : 1000
dmc_stats_nstep : 10000
EOF

casino > dmc.out
```

### Using runqmc:
```bash
# Automated workflow
runqmc --define "nproc=16" --define "system=material"
```

## Advanced Features

### Geminal Wavefunctions:
- Pfaffian determinants
- Pairing functions
- BCS-like states
- Superconductors
- Advanced correlation

### Backflow:
- Electron backflow
- Enhanced trial functions
- Lower variance
- Better nodes
- Systematic improvement

### Released-Node DMC:
- Beyond fixed-node
- Improved accuracy
- Sign problem mitigation
- Advanced method

### Excited States:
- Multiple states
- Promotion methods
- Optical properties
- Excitation energies

## Performance Characteristics
- **Speed**: Efficient, MPI parallel
- **Accuracy**: Benchmark quality
- **System size**: Moderate to large
- **Purpose**: Research and production
- **Typical**: HPC calculations

## Computational Cost
- DMC expensive
- Advanced features costly
- HPC recommended
- Production capable
- Benchmark quality justifies cost

## Limitations & Known Constraints
- **Academic license**: Registration required
- **Not fully open**: Source available to academics
- **Learning curve**: Steep
- **Fixed-node**: DMC nodal approximation
- **Computational cost**: Expensive
- **Documentation**: Extensive but complex

## Comparison with Other QMC Codes
- **vs QMCPACK**: CASINO feature-rich, QMCPACK HPC-optimized
- **vs TurboRVB**: CASINO general, TurboRVB specialized
- **Unique strength**: Mature features, geminals, backflow, methodological depth, community standard

## Application Areas

### Electronic Structure:
- Molecules
- Crystals
- Surfaces
- 2D materials
- Nanostructures

### Benchmark Calculations:
- High-accuracy reference
- Method validation
- DFT comparison
- Chemical accuracy

### Method Development:
- Advanced wavefunctions
- New QMC methods
- Algorithm research
- Methodological studies

### Materials Science:
- Energetics
- Structures
- Excited states
- Properties

## Best Practices

### Trial Wavefunctions:
- Quality DFT starting point
- Geminals for pairing
- Backflow optimization
- Multi-reference when needed

### DMC:
- Timestep extrapolation
- Population control
- Finite-size corrections
- Careful error analysis

### Feature Usage:
- Understand advanced features
- Read documentation thoroughly
- Start simple
- Validate results

## Community and Support
- Academic license (free for academics)
- Large user community
- Active mailing list
- Documentation
- Workshops
- Publications

## Educational Resources
- Comprehensive manual
- Tutorials
- Example inputs
- QMC school materials
- Publication list
- User forum

## Development
- Cambridge origin
- International collaboration
- Active development
- Regular updates
- Community contributions
- Method innovation

## Research Impact
CASINO has enabled numerous high-impact QMC studies and method development over decades, serving as a community standard for advanced QMC calculations with thousands of publications.

## Verification & Sources
**Primary sources**:
1. Homepage: https://vallico.net/casino/
2. Documentation
3. Publications: J. Chem. Phys. 152, 154106 (2020)

**Secondary sources**:
1. QMC literature
2. User publications
3. Method development papers

**Confidence**: CONFIRMED - Established QMC code

**Verification status**: ✅ CONFIRMED
- Website: ACTIVE
- License: Academic (free for academics)
- **Category**: Academic QMC code
- Status: Actively developed
- Community: Large, international
- Specialized strength: Feature-rich quantum Monte Carlo, VMC/DMC methods, geminal wavefunctions, backflow transformations, released-node DMC, decades of development, community standard, benchmark calculations, advanced trial functions, comprehensive methodology, mature production code
