# Gollum

## Official Resources
- Source Repository: https://github.com/gollumcode/gollum2
- Documentation: Included in repository
- License: Open source

## Overview
**Gollum** is a next-generation quantum transport simulation tool for computing transport properties of nanoscale devices using the non-equilibrium Green's function (NEGF) method. It works with tight-binding Hamiltonians and can compute conductance, current-voltage characteristics, and thermoelectric properties.

**Scientific domain**: Quantum transport, NEGF, molecular electronics  
**Target user community**: Researchers simulating quantum transport in molecular junctions and nanoscale devices

## Theoretical Methods
- Non-equilibrium Green's function (NEGF) method
- Landauer-Büttiker formalism
- Tight-binding Hamiltonians
- Density functional theory (DFT) input
- Self-energy calculation
- Transmission function calculation

## Capabilities (CRITICAL)
- Quantum conductance calculation
- Current-voltage (I-V) characteristics
- Transmission function
- Thermoelectric coefficients
- Tight-binding model transport
- DFT Hamiltonian input
- Multi-terminal transport
- Spin-dependent transport

**Sources**: GitHub repository, J. Chem. Phys.

## Key Strengths

### NEGF Transport:
- Full Green's function calculation
- Self-energy for leads
- Bias-dependent transport
- Multi-terminal support

### Versatile Hamiltonian:
- Tight-binding models
- DFT-derived Hamiltonians
- Custom Hamiltonians
- Parameter exploration

### Thermoelectric:
- Seebeck coefficient
- Peltier coefficient
- Thermal conductance
- ZT figure of merit

## Inputs & Outputs
- **Input formats**:
  - Hamiltonian matrix files
  - Lead self-energy parameters
  - Transport configuration
  
- **Output data types**:
  - Transmission vs energy
  - I-V characteristics
  - Thermoelectric coefficients
  - Local density of states

## Interfaces & Ecosystem
- **DFT codes**: Hamiltonian extraction
- **Python**: Scripting
- **NumPy**: Numerical computation

## Performance Characteristics
- **Speed**: Fast for TB models
- **Accuracy**: Depends on Hamiltonian quality
- **System size**: Thousands of orbitals
- **Memory**: Moderate

## Computational Cost
- **Transmission**: Seconds to minutes
- **I-V curve**: Minutes
- **Typical**: Efficient

## Limitations & Known Constraints
- **Tight-binding**: Quality depends on Hamiltonian
- **No self-consistent NEGF**: No Poisson-NEGF
- **Limited documentation**: Research code
- **Small community**: Research group code

## Comparison with Other Codes
- **vs Transiesta**: Gollum is TB/NEGF, Transiesta is DFT+NEGF
- **vs Kwant**: Gollum focuses on molecular junctions, Kwant is general TB
- **vs Nanodcal**: Gollum is open source, Nanodcal is commercial
- **Unique strength**: Next-generation quantum transport for molecular junctions, thermoelectric properties

## Application Areas

### Molecular Electronics:
- Molecular junction conductance
- Single-molecule transport
- Break junction simulations
- Switching behavior

### Thermoelectrics:
- Molecular thermoelectrics
- Seebeck coefficient prediction
- ZT optimization
- Energy harvesting

### Nanoscale Devices:
- Quantum dot transport
- Nanowire conductance
- 2D material junctions
- Spin-dependent transport

## Best Practices

### Hamiltonian Quality:
- Use well-converged DFT Hamiltonians
- Test convergence with basis size
- Validate against known systems
- Consider spin-orbit coupling

### Transport Calculation:
- Use sufficient energy grid
- Include sufficient lead layers
- Check convergence of self-energies
- Compare with Landauer limit

## Community and Support
- Open source on GitHub
- Research code
- Limited documentation
- Related publications available

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/gollumcode/gollum2
2. Related publications from Lancaster University

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Research code
- Specialized strength: Quantum transport for molecular junctions, NEGF, thermoelectric properties
