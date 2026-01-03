# TurboRVB (Resonating Valence Bond Quantum Monte Carlo)

## Official Resources
- Homepage: https://turborvb.qe-forge.org/
- Documentation: https://turborvb.qe-forge.org/documentation/
- Source Repository: https://github.com/sissaschool/turborvb
- License: GNU General Public License v3.0

## Overview
TurboRVB is a high-performance quantum Monte Carlo package developed at SISSA (International School for Advanced Studies, Trieste) with emphasis on strongly correlated systems, superconductors, and resonating valence bond (RVB) physics. The code implements advanced trial wavefunctions including Jastrow-geminal-Slater forms, AGP (antisymmetrized geminal power), and pairing functions optimized for studying correlation effects, superconductivity, and quantum phase transitions. TurboRVB is GPU-accelerated and designed for large-scale calculations.

**Scientific domain**: Quantum Monte Carlo, strongly correlated systems, superconductivity  
**Target user community**: Correlated materials researchers, superconductivity studies, QMC specialists

## Theoretical Methods
- Variational Monte Carlo (VMC)
- Lattice regularized diffusion Monte Carlo (LRDMC)
- Jastrow-geminal-Slater wavefunctions
- AGP (antisymmetrized geminal power)
- Pairing wavefunctions
- Pfaffian determinants
- Resonating valence bond states
- BCS-like wavefunctions
- GPU-accelerated algorithms

## Capabilities (CRITICAL)
**Category**: Open-source QMC code (specialized)
- VMC and LRDMC methods
- Advanced pairing wavefunctions
- Geminal functions
- AGP ansatz
- Strongly correlated systems
- Superconductors
- Hubbard models
- Periodic systems
- GPU acceleration (CUDA)
- Wavefunction optimization
- Energy and forces
- Excited states
- Production quality

**Sources**: Official website, GitHub, publications

## Key Strengths

### Pairing Wavefunctions:
- AGP ansatz
- Geminal functions
- BCS-like states
- Pfaffian forms
- Superconductivity-optimized

### Strongly Correlated:
- Designed for correlations
- Hubbard models
- Quantum magnetism
- Mott physics
- RVB states

### GPU Performance:
- CUDA acceleration
- High throughput
- Large-scale systems
- Optimized kernels
- Modern HPC

### SISSA Development:
- Expert group
- Research-driven
- Active development
- Superconductivity focus
- Method innovation

## Inputs & Outputs
- **Input formats**:
  - TurboRVB input files
  - DFT trial wavefunctions
  - Lattice model definitions
  - Wavefunction parameters
  
- **Output data types**:
  - Total energies
  - Pairing correlations
  - Order parameters
  - Forces
  - Observables
  - Wavefunction data

## Interfaces & Ecosystem

### DFT Integration:
- Quantum ESPRESSO
- Trial wavefunction input
- Real-space projections

### GPU Computing:
- CUDA support
- Mixed precision
- Performance optimization
- Large-scale capability

## Workflow and Usage

### Installation:
```bash
# Clone repository
git clone https://github.com/sissaschool/turborvb.git
cd turborvb
# Configure and build (with GPU)
./configure --enable-parallel --enable-gpu
make
```

### Lattice Model:
```bash
# Hubbard model example
# Define lattice and parameters
# Setup AGP wavefunction
# Run VMC optimization
```

### Wavefunction Optimization:
```bash
# Optimize pairing function
turborvb-optimize.x < input.d
```

### Energy Calculation:
```bash
# VMC or LRDMC
turborvb.x < input.d
```

## Advanced Features

### AGP Wavefunction:
- Antisymmetrized geminal power
- Pairing ansatz
- Superconducting correlations
- Optimized for BCS-like states
- Pfaffian evaluation

### Geminal Functions:
- General pairing
- Flexible correlations
- Beyond Slater determinants
- Advanced trial functions

### Strongly Correlated Models:
- Hubbard model
- t-J model
- Extended Hubbard
- Quantum magnetism
- Lattice systems

### LRDMC:
- Lattice regularized DMC
- Alternative to standard DMC
- Specific advantages
- Production quality

## Performance Characteristics
- **Speed**: GPU-accelerated, fast
- **Accuracy**: High quality
- **System size**: Large systems (GPU)
- **Purpose**: Correlated systems, superconductivity
- **Typical**: GPU workstations to HPC

## Computational Cost
- GPU acceleration crucial
- Efficient for large systems
- Wavefunction optimization expensive
- Production capable
- HPC-suitable

## Limitations & Known Constraints
- **Specialized focus**: Pairing/correlations
- **GPU recommended**: Best performance
- **Learning curve**: Advanced wavefunctions
- **Documentation**: Growing
- **Community**: Smaller than QMCPACK/CASINO
- **Trial functions**: Requires expertise

## Comparison with Other QMC Codes
- **vs QMCPACK**: TurboRVB pairing-specialized, QMCPACK general
- **vs CASINO**: TurboRVB GPU-focused, CASINO feature-rich
- **Unique strength**: AGP/geminal wavefunctions, superconductivity, GPU performance, strongly correlated focus, RVB physics

## Application Areas

### Superconductivity:
- High-Tc materials
- Pairing mechanisms
- BCS vs exotic pairing
- Order parameters
- Phase transitions

### Strongly Correlated:
- Hubbard model
- Mott insulators
- Quantum magnetism
- Correlation effects
- Phase diagrams

### Quantum Materials:
- Cuprates
- Pnictides
- Correlated electrons
- Quantum criticality
- Exotic phases

## Best Practices

### Wavefunction Choice:
- AGP for pairing systems
- Geminals for correlations
- Start simple, add complexity
- Systematic optimization

### GPU Usage:
- CUDA-enabled GPUs
- Mixed precision
- Performance tuning
- Resource optimization

### Optimization:
- Careful wavefunction setup
- Parameter optimization
- Convergence testing
- Physical validation

## Community and Support
- Open-source (GPL v3)
- SISSA development
- GitHub repository
- Research community
- Growing user base
- Scientific publications

## Educational Resources
- Official documentation
- GitHub examples
- Scientific papers
- SISSA workshops
- User contributions

## Development
- SISSA (Trieste, Italy)
- Active research group
- GPU focus
- Method development
- Superconductivity expertise
- Regular updates

## Research Impact
TurboRVB enables advanced QMC studies of strongly correlated systems and superconductors, particularly valuable for exploring pairing mechanisms and exotic correlation effects beyond standard trial wavefunctions.

## Verification & Sources
**Primary sources**:
1. Homepage: https://turborvb.qe-forge.org/
2. GitHub: https://github.com/sissaschool/turborvb
3. Publications: J. Chem. Phys. 152, 204121 (2020)

**Secondary sources**:
1. QMC literature
2. Superconductivity papers
3. User publications

**Confidence**: CONFIRMED - Specialized QMC code

**Verification status**: ✅ CONFIRMED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: GPL v3 (open-source)
- **Category**: Open-source QMC code
- Status: Actively developed
- Institution: SISSA (Trieste)
- Specialized strength: Advanced pairing wavefunctions (AGP/geminal), superconductivity studies, strongly correlated systems, GPU acceleration, resonating valence bond physics, Pfaffian forms, LRDMC method, SISSA development, specialized for correlation and pairing physics
