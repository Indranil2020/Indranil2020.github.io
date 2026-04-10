# Jiezi

## Official Resources
- Source Repository: https://github.com/Jiezi-negf/Jiezi
- Documentation: Included in repository
- License: Open source

## Overview
**Jiezi** is an open-source Python software for simulating quantum transport of nanoscale devices. It solves the Schrödinger and Poisson equations self-consistently using the non-equilibrium Green's function (NEGF) method with finite element discretization, enabling atomistic-level device simulation.

**Scientific domain**: Self-consistent NEGF quantum transport, device simulation  
**Target user community**: Researchers simulating quantum transport in nanoscale electronic devices with self-consistent electrostatics

## Theoretical Methods
- Non-equilibrium Green's function (NEGF)
- Self-consistent Schrödinger-Poisson solver
- Finite element method (FEM)
- Tight-binding Hamiltonians
- Density matrix calculation
- Open boundary conditions
- Scattering self-energies

## Capabilities (CRITICAL)
- Self-consistent NEGF-Poisson simulation
- Quantum transport in nanoscale devices
- I-V characteristics with self-consistent potential
- Transmission function
- Density of states
- Charge density distribution
- Electrostatic potential
- FET and nanowire device simulation

**Sources**: GitHub repository

## Key Strengths

### Self-Consistent NEGF:
- Schrödinger-Poisson coupling
- Realistic device electrostatics
- Gate voltage effects
- Charge self-consistency

### Finite Element Method:
- Flexible geometry
- Non-uniform meshing
- Accurate potential
- Complex device shapes

### Python Framework:
- Easy to use and modify
- Jupyter notebook compatible
- Extensible architecture
- Clear code structure

## Inputs & Outputs
- **Input formats**:
  - Python configuration scripts
  - Device geometry parameters
  - Material parameters
  
- **Output data types**:
  - I-V characteristics
  - Transmission spectra
  - Charge density
  - Electrostatic potential
  - Density of states

## Interfaces & Ecosystem
- **Python**: Primary language
- **NumPy/SciPy**: Numerical computation
- **Matplotlib**: Visualization
- **FEM**: Finite element discretization

## Performance Characteristics
- **Speed**: Moderate (self-consistent iteration)
- **Accuracy**: Good with converged Hamiltonian
- **System size**: Thousands of atoms
- **Memory**: Moderate to high

## Computational Cost
- **Single bias point**: Minutes to hours
- **Full I-V**: Hours to days
- **Typical**: Moderate to expensive

## Limitations & Known Constraints
- **Tight-binding only**: No DFT integration
- **Python speed**: Slower than compiled codes
- **Documentation**: Limited
- **Small community**: Research code

## Comparison with Other Codes
- **vs Transiesta**: Jiezi is Python/TB, Transiesta is DFT+NEGF
- **vs Nanodcal**: Jiezi is open source, Nanodcal is commercial LCAO
- **vs Smeagol**: Jiezi is standalone Python, Smeagol is SIESTA-based
- **Unique strength**: Self-consistent NEGF-Poisson with FEM, Python framework for device simulation

## Application Areas

### Nanoscale FETs:
- Nanowire FET simulation
- Gate voltage effects
- Subthreshold swing
- ON/OFF current

### Molecular Devices:
- Molecular junction I-V
- Gate-controlled transport
- Single-molecule transistors
- Switching behavior

### 2D Material Devices:
- Graphene nanoribbon FETs
- MoS2 transistor simulation
- Heterostructure devices
- Contact resistance

## Best Practices

### Self-Consistent Convergence:
- Monitor charge convergence
- Use appropriate mixing
- Start from zero-bias solution
- Check Poisson convergence

### Device Setup:
- Use realistic geometry
- Include sufficient lead length
- Appropriate boundary conditions
- Validate against analytical models

## Community and Support
- Open source on GitHub
- Research code
- Limited documentation
- Example calculations provided

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/Jiezi-negf/Jiezi

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Research code
- Specialized strength: Self-consistent NEGF-Poisson with FEM, Python device simulation
