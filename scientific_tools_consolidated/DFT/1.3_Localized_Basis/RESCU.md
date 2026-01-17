# RESCU (Real-space Electronic Structure CalcUlator)

## Official Resources
- Homepage: https://www.nanoacademic.com/rescu
- Documentation: https://www.nanoacademic.com/rescu/documentation
- Developer: NanoAcademic Technologies / McGill University
- License: Commercial

## Overview
RESCU is a Kohn-Sham DFT solver that combines multiple basis set approaches (atomic orbitals, plane waves, real-space grids) within a single framework. Developed primarily in MATLAB with C extensions, it is designed for large-scale simulations of materials containing thousands to tens of thousands of atoms with modest computational resources.

**Scientific domain**: Semiconductors, 2D materials, nanoelectronics, photovoltaics  
**Target user community**: Materials scientists and device engineers requiring large-scale DFT for realistic nanostructures

## Theoretical Methods
- Density Functional Theory (DFT)
- Numerical Atomic Orbitals (NAOs)
- Plane-wave expansion
- Real-space grid discretization
- LDA, GGA, meta-GGA functionals
- Hartree-Fock exchange
- Hybrid functionals (HSE, PBE0)
- DFT+U for correlated systems
- Modified Becke-Johnson (mBJ)

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Large-scale calculations (1000s-10000s atoms)
- Band structure and DOS
- Projected DOS (PDOS)
- Optical properties
- Multiple boundary conditions
- Geometry optimization
- Electronic transport (NEGFwith QuantumATK)
- k-point sampling
- Spin-polarized calculations
- Parallel execution

**Sources**: NanoAcademic Technologies, arXiv publications

## Key Strengths

### Hybrid Basis Approach:
- Combines NAO, plane-wave, real-space
- Single framework flexibility
- Method benchmarking capability
- Optimal for different systems

### Large-Scale Capability:
- Tens of thousands of atoms
- Modest computational resources
- Linear scaling techniques
- Efficient memory management

### Delayed Cubic Scaling:
- O(N³) onset at larger sizes
- NAO initial subspace
- Occupied-state focus
- Efficient for device simulations

### Functional Diversity:
- Semi-local (LDA, GGA, mGGA)
- Hybrid (HSE, PBE0)
- DFT+U
- mBJ for band gaps

## Inputs & Outputs
- **Input formats**:
  - Structure files
  - MATLAB interface
  - Parameter specifications
  - Pseudopotential files
  
- **Output data types**:
  - Total energies
  - Band structure
  - DOS/PDOS
  - Charge densities
  - Wave functions
  - Optical spectra

## Interfaces & Ecosystem
- **QuantumATK integration**:
  - Transport calculations
  - Device simulations
  - Graphical interface
  
- **Analysis tools**:
  - MATLAB post-processing
  - Visualization scripts
  - Property extraction

## Advanced Features

### Multi-Method Capability:
- Switch between NAO, PW, real-space
- Systematic comparison
- Method validation
- Research flexibility

### Device-Scale DFT:
- Thousands of atoms routine
- Realistic nanostructures
- Heterostructure modeling
- Interface calculations

### Hybrid Functionals:
- Efficient implementation
- Large system hybrids
- Accurate band gaps
- Electronic properties

### Transport Coupling:
- Interface to NEGF methods
- Device simulations
- Current calculations
- Quantum transport

## Performance Characteristics
- **Speed**: Efficient for large systems
- **Accuracy**: Standard DFT accuracy
- **System size**: Up to tens of thousands atoms
- **Memory**: Optimized management
- **Parallelization**: Multi-core and distributed

## Computational Cost
- **Large systems**: Efficient delayed scaling
- **Hybrid DFT**: Feasible for large systems
- **Typical**: Workstation to cluster
- **Memory**: Careful management for size

## Limitations & Known Constraints
- **Commercial license**: Not freely available
- **MATLAB core**: Requires MATLAB license
- **Specialization**: Materials focus
- **Community**: Smaller than major codes
- **Documentation**: Commercial-level

## Comparison with Other Codes
- **vs VASP**: RESCU hybrid basis vs VASP plane-wave
- **vs SIESTA**: Both NAO-capable, different architectures
- **vs QuantumATK**: RESCU integrates, different focuses
- **Unique strength**: Multi-basis hybrid, large-scale device DFT

## Application Areas

### Nanoelectronics:
- Transistor channels
- 2D material devices
- Heterostructure electronics
- Contact interfaces

### Photovoltaics:
- Solar cell materials
- Interfaces and junctions
- Optical absorption
- Defects in absorbers

### 2D Materials:
- Graphene nanostructures
- TMD devices
- Heterostructure stacking
- Edge effects

### Semiconductor Devices:
- Realistic device regions
- Source-drain channels
- Gate interfaces
- Doping profiles

## Best Practices

### Basis Selection:
- NAO for initial efficiency
- Plane-wave for validation
- Match to system type

### Large Systems:
- Use NAO mode primarily
- Optimize k-points
- Monitor memory usage

### Hybrid Functionals:
- Test on smaller systems first
- Balance accuracy and cost

## Community and Support
- NanoAcademic Technologies
- Commercial support
- Training courses
- Academic publications
- McGill University development

## Verification & Sources
**Primary sources**:
1. NanoAcademic: https://www.nanoacademic.com/rescu
2. arXiv: RESCU methodology papers
3. M. Côté group publications (McGill)

**Confidence**: VERIFIED - Commercial product, published methodology

**Verification status**: ✅ VERIFIED
- Source code: Commercial
- Academic use: Publications with RESCU
- Documentation: Commercial quality
- Active development: Commercial updates
- Specialty: Large-scale DFT, multi-basis hybrid, device simulations
