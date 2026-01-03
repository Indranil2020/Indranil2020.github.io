# NECI (N-Electron Configuration Interaction)

## Official Resources
- Homepage: https://github.com/ghb24/NECI_STABLE
- Documentation: GitHub repository and manual
- Source Repository: https://github.com/ghb24/NECI_STABLE
- License: GNU General Public License v3.0

## Overview
NECI is a state-of-the-art implementation of Full Configuration Interaction Quantum Monte Carlo (FCIQMC), a stochastic method for solving the electronic Schrödinger equation in a systematically improvable way. Developed primarily at the University of Cambridge, NECI provides numerically exact solutions to the many-electron problem by stochastically sampling the full CI space. The method bridges quantum chemistry and QMC, offering chemical accuracy for strongly correlated systems.

**Scientific domain**: FCIQMC, quantum chemistry, strongly correlated electrons  
**Target user community**: Quantum chemists, strongly correlated systems, method developers

## Theoretical Methods
- Full Configuration Interaction QMC (FCIQMC)
- Initiator approximation
- Spawning algorithm
- Stochastic CI sampling
- Coupled Cluster Monte Carlo (CCMC)
- Density Matrix QMC (DMQMC)
- Finite-temperature extensions
- Exact (in CI space limit)

## Capabilities (CRITICAL)
**Category**: Open-source FCIQMC code
- FCIQMC method
- Initiator FCIQMC
- CCMC variants
- DMQMC (finite-T)
- Molecules and lattices
- Strongly correlated systems
- Exact in complete space
- Systematically improvable
- Chemical accuracy
- Ground and excited states
- MPI parallelization
- Production quality

**Sources**: GitHub repository, Cambridge group publications

## Key Strengths

### FCIQMC Method:
- Stochastic full CI
- Numerically exact (converged)
- Strong correlation capable
- Systematically improvable
- Sign coherence

### Initiator Approximation:
- Controlled approximation
- Larger systems
- Sign problem mitigation
- Systematic convergence
- Production viable

### Versatility:
- FCIQMC, CCMC, DMQMC variants
- Molecules and lattices
- Multiple methods
- Research flexibility
- Method development

### Cambridge Development:
- Leading research group
- Method innovation
- Active development
- Strong theoretical foundation
- Publication quality

## Inputs & Outputs
- **Input formats**:
  - NECI input files
  - Integrals from quantum chemistry codes
  - FCIDUMP format
  - Model Hamiltonians
  
- **Output data types**:
  - Ground state energies
  - Excited states
  - Correlation energies
  - Observables
  - Sampling statistics

## Interfaces & Ecosystem

### Quantum Chemistry:
- FCIDUMP integrals
- Molpro
- PySCF
- GAMESS
- Standard formats

### Lattice Models:
- Hubbard model
- Custom Hamiltonians
- Research applications

## Workflow and Usage

### Installation:
```bash
# Clone repository
git clone https://github.com/ghb24/NECI_STABLE.git
cd NECI_STABLE
mkdir build && cd build
cmake ..
make -j8
```

### Input File (FCIQMC.inp):
```
FCIQMC
  METHODS
    method vertex
  ENDINIT
  
  CALC
    electrons 10
    spin-restrict
    totalwalkers 1e6
    startsinglepart 100
  ENDCALC
  
  LOGGING
    popsfile-format HDF5
  ENDLOG
END
```

### Run FCIQMC:
```bash
# MPI parallel
mpirun -n 16 neci FCIQMC.inp > output.out
```

### Analysis:
```bash
# Extract energy
grep "Shift" output.out | tail -100
```

## Advanced Features

### Initiator Approximation:
- i-FCIQMC
- Controlled bias
- Sign coherence
- Larger systems
- Systematic convergence

### CCMC:
- Coupled Cluster Monte Carlo
- Truncated CC
- Alternative approach
- Specific advantages

### DMQMC:
- Density Matrix QMC
- Finite temperature
- Thermal properties
- Imaginary-time evolution

### Excited States:
- Multiple states
- State-specific
- Systematic approach
- Excitation energies

## Performance Characteristics
- **Speed**: MPI-parallel
- **Accuracy**: Numerically exact (converged)
- **System size**: Moderate (CI space limited)
- **Purpose**: Strongly correlated, benchmarks
- **Typical**: HPC calculations

## Computational Cost
- Walker-number dependent
- CI space scaling
- Expensive for large systems
- Exact results justify cost
- HPC recommended

## Limitations & Known Constraints
- **CI space**: Limited by basis
- **Computational cost**: Expensive
- **Sign problem**: Fermion sign issue
- **System size**: Moderate
- **Learning curve**: FCIQMC expertise
- **HPC required**: Production calculations

## Comparison with Other Methods
- **vs Traditional FCI**: NECI stochastic, tractable for larger
- **vs DMRG**: NECI different approach, complementary
- **vs Coupled Cluster**: NECI exact, CC approximate
- **Unique strength**: FCIQMC method, systematically exact, strong correlation, Cambridge development

## Application Areas

### Strongly Correlated:
- Molecules
- Transition metals
- f-electron systems
- Strong correlation
- Multi-reference

### Benchmark Calculations:
- Exact results
- Method validation
- Chemical accuracy
- Reference energies
- Correlation energies

### Quantum Chemistry:
- Ground states
- Excited states
- Reaction energies
- Spectroscopy
- Electronic structure

### Method Development:
- FCIQMC research
- Algorithm development
- Stochastic methods
- QMC innovations

## Best Practices

### Walker Number:
- Sufficient population
- Convergence testing
- Statistical analysis
- Plateau region
- Error estimation

### Initiator:
- Appropriate threshold
- Systematic convergence
- Balance accuracy/cost
- Production settings

### Basis Sets:
- Quality basis
- Convergence studies
- CI space considerations
- System-appropriate

## Community and Support
- Open-source (GPL v3)
- University of Cambridge
- GitHub repository
- Active development
- Research community
- Publications
- Method innovation

## Educational Resources
- GitHub documentation
- FCIQMC papers
- User manual
- Example inputs
- Cambridge group publications
- QMC schools

## Development
- University of Cambridge
- Alavi group
- Active research
- Method development
- Community contributions
- Regular updates

## Research Impact
NECI and FCIQMC have revolutionized exact solutions for strongly correlated molecules, enabling benchmark-quality calculations that were previously impossible, with hundreds of publications and major impact in quantum chemistry.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ghb24/NECI_STABLE
2. Cambridge group
3. Publications: J. Chem. Phys. 132, 041103 (2010)

**Secondary sources**:
1. FCIQMC literature
2. Quantum chemistry papers
3. User publications

**Confidence**: VERIFIED - Leading FCIQMC code

**Verification status**: ✅ VERIFIED
- GitHub: ACCESSIBLE
- License: GPL v3 (open-source)
- **Category**: Open-source FCIQMC code
- Status: Actively developed
- Institution: University of Cambridge
- Specialized strength: Full Configuration Interaction Quantum Monte Carlo, stochastic CI sampling, numerically exact for strongly correlated systems, initiator approximation, CCMC/DMQMC variants, benchmark quality, Cambridge development, chemical accuracy, systematic improvability, production quality
