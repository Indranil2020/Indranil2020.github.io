# QuaTrEx24

## Official Resources
- Homepage: Research code (referenced in arXiv:2508.19138)
- Documentation: "Ab-initio Quantum Transport with the GW Approximation..." (arXiv:2508.19138)
- Source Repository: Contact authors (ETH Zurich / UT Austin / related benchmarks)
- License: Research / Open Source availability pending

## Overview
QuaTrEx24 is a cutting-edge research code implementing an optimized Non-Equilibrium Green's Function (NEGF) method combined with the GW approximation (NEGF+GW). It is designed for sustained exascale performance on leadership-class supercomputers (e.g., Alps, Frontier), demonstrating calculations on realistic device structures with over 42,000 atoms (and scaling up to ~85,000 atoms).

**Scientific domain**: Quantum transport, nanoscale devices, large-scale GW, NEGF, Exascale computing  
**Target user community**: Device physicists, HPC specialists, transistor modeling researchers

## Theoretical Methods
- NEGF formalism
- GW approximation
- Optimized sparse algebra
- Large-scale massively parallel solvers
- Finite-difference / localized basis
- Quantum transport properties
- Electron-electron correlation in devices

## Capabilities (CRITICAL)
- NEGF+GW transport calculations
- Massive scale systems (42,240+ atoms demonstrated)
- Nanowire and transistor simulations
- Sustained exascale performance
- Electronic current/transport with correlations
- Quasiparticle corrections in operating devices
- State-of-the-art scaling (2025 benchmark)

**Sources**: arXiv:2508.19138 "Ab-initio Quantum Transport with the GW Approximation, 42,240 Atoms..."

## Key Strengths

### Exascale Scaling:
- Demonstrated on >40,000 atoms
- Scales to full machine partition (Frontier/Alps)
- Optimized sparse algorithms
- Breakthrough capacity for GW transport

### Transport Focus:
- Direct NEGF integration
- Transport + Correlation
- Beyond standard band structure
- Device characteristics

### Modern Optimization:
- 2024 algorithmic advances
- Efficient memory usage
- Advanced solver techniques

## Inputs & Outputs
- **Input formats**:
  - Device geometry
  - Hamiltonian data
  - Transport parameters
  
- **Output data types**:
  - Transmission functions
  - Current-voltage characteristics
  - Spectral functions
  - GW renormalized levels

## Interfaces & Ecosystem
- **Integration**: Likely standalone or specific transport framework
- **Basis**: Localized/Tight-binding typically
- **Status**: Research grade

## Performance Characteristics
- **Speed**: Highly optimized for size
- **Accuracy**: GW level in transport
- **System size**: Massive (10k atoms)
- **Parallelization**: HPC ready

## Computational Cost
- **High absolute cost**: Large systems
- **Efficiency**: High per-atom efficiency
- **Feasibility**: Enables previously impossible sizes

## Limitations & Known Constraints
- **Availability**: Research code status
- **Documentation**: Limited to publications
- **Generalization**: focused on devices

## Comparison with Other Codes
- **vs BerkeleyGW/Yambo**: QuaTrEx24 focuses on *transport* (NEGF) and massive scaling
- **vs OpenAtom**: Different methodology (NEGF vs Real-space)
- **Unique strength**: NEGF+GW at 10,000 atom scale

## Application Areas

### Nanoelectronics:
- Transistor modeling
- Nanowire devices
- Quantum transport
- Device scaling limits

### Emerging Devices:
- 2D material devices
- Molecular electronics
- Correlated transport

## Community and Support
- Research group maintained
- Publication-based support
- Cutting-edge academic research

## Verification & Sources
**Primary sources**:
1. arXiv preprint (Aug 2025/2024 activity)
2. Simulation mentions in recent transport literature

**Confidence**: MEDIUM (Research Code)
- Existence: CONFIRMED via literature
- Availability: Research/Contact authors
