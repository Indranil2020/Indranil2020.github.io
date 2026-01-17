# Dice

## Official Resources
- Homepage: https://sanshar.github.io/Dice/
- Documentation: https://sanshar.github.io/Dice/
- Source Repository: https://github.com/sanshar/Dice
- License: Apache License 2.0

## Overview
Dice is a quantum chemistry software that implements the Semistochastic Heat-Bath Configuration Interaction (SHCI) algorithm. It provides a highly efficient approach to near-Full Configuration Interaction (FCI) calculations for systems with large active spaces. Developed by Sandeep Sharma's group at the University of Colorado Boulder, Dice can handle active spaces of 30-100 orbitals, far beyond traditional CASCI/CASSCF limits.

**Scientific domain**: Strong correlation, selected configuration interaction, near-FCI methods  
**Target user community**: Researchers studying strongly correlated systems requiring large active space treatments

## Theoretical Methods
- Semistochastic Heat-Bath Configuration Interaction (SHCI)
- Selected Configuration Interaction
- Perturbative corrections (PT2)
- Variational CI optimization
- Stochastic sampling methods
- Determinant selection algorithms
- Multireference treatment

## Capabilities (CRITICAL)
- Near-FCI accuracy for large active spaces
- Active spaces of 30-100 orbitals
- Heat-bath selection algorithm
- Semistochastic PT2 corrections
- Efficient determinant generation
- Variational optimization
- Memory-efficient algorithms
- Interface with PySCF and other codes
- Ground and excited states
- Spin-pure wavefunctions

## Key Strengths

### Large Active Spaces:
- 30-100 orbital active spaces
- Beyond CASCI/CASSCF limits
- Near-FCI accuracy
- Polynomial scaling for selected CI

### SHCI Algorithm:
- Efficient determinant selection
- Heat-bath criterion
- Stochastic PT2
- Low memory requirements
- Unbiased selection

### Strong Correlation:
- Multi-reference problems
- Transition metal clusters
- Bond breaking
- Biradical systems
- Open-shell singlets

### Integration:
- PySCF interface
- Integral file input
- Molpro/Gaussian compatibility
- Script-driven workflow

## Inputs & Outputs
- **Input formats**:
  - FCIDUMP files
  - PySCF interface
  - Custom input files
  - Integral files
  
- **Output data types**:
  - Variational energies
  - PT2-corrected energies
  - CI coefficients
  - Excitation energies
  - Natural occupation numbers

## Interfaces & Ecosystem
- **PySCF integration**: Primary interface
- **Molpro**: FCIDUMP generation
- **Gaussian**: Integral conversion
- **SHCI-PDFT**: Interface for PDFT
- **QMC codes**: Trial wavefunction generation

## Advanced Features

### Heat-Bath Selection:
- Importance sampling
- Efficient screening
- Deterministic core
- Stochastic tail
- Convergence control

### Excited States:
- State-averaged calculations
- Multiple roots
- Transition moments
- Spin state selection

### Perturbative Corrections:
- Semistochastic PT2
- Low-memory algorithms
- Parallelizable
- Accurate extrapolation

### Analysis Tools:
- Natural orbital occupations
- Correlation diagnostics
- Determinant analysis
- Spin populations

## Performance Characteristics
- **Speed**: Highly efficient for selected CI
- **Accuracy**: Near-FCI quality
- **System size**: Large active spaces (30-100 orbitals)
- **Memory**: Efficient algorithms
- **Parallelization**: MPI parallelized

## Computational Cost
- **Variational SHCI**: Polynomial in selected determinants
- **PT2 correction**: Additional overhead
- **Scaling**: Better than FCI for large systems
- **Typical**: Hours to days for complex systems

## Limitations & Known Constraints
- **Dynamic correlation**: PT2 approximate
- **Total electrons**: Still limited by active space
- **Properties**: Energy-focused
- **Learning curve**: Specialized method
- **Post-processing**: Limited analysis tools

## Comparison with Other Codes
- **vs CASSCF**: Dice handles larger active spaces
- **vs DMRG**: Different correlation approach
- **vs FCIQMC**: Deterministic core, different sampling
- **vs Quantum Package CIPSI**: Different selection criteria
- **Unique strength**: Large active spaces with near-FCI accuracy

## Application Areas

### Transition Metal Chemistry:
- Metal clusters
- Spin state energetics
- Multicenter bonding
- Catalytic intermediates

### Strong Correlation:
- Bond dissociation
- Biradicals
- Polyradicals
- Open-shell singlets

### Excited States:
- Electronic spectra
- State orderings
- Avoided crossings
- Conical intersections

### Benchmarking:
- Reference calculations
- Method validation
- Correlation analysis
- Near-exact energies

## Best Practices

### Active Space Selection:
- Include correlating orbitals
- Balanced occupied/virtual
- Test convergence with size
- Natural orbital analysis

### Convergence:
- Monitor variational energy
- Check PT2 contribution
- Convergence with determinants
- Multiple runs for stochastic error

### Calculation Setup:
- Start from good orbitals
- Use FCIDUMP format
- Appropriate memory settings
- Parallel execution

## Community and Support
- Open source Apache 2.0
- Active GitHub development
- Sharma group (CU Boulder)
- Academic publications
- Growing user community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/sanshar/Dice
2. Documentation: https://sanshar.github.io/Dice/
3. Holmes et al., J. Chem. Theory Comput. 12, 3674 (2016) - SHCI paper
4. Sharma et al., J. Chem. Phys. 149, 214110 (2018)

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, Apache 2.0)
- Documentation: Available
- Active development: Yes
- Widely cited: High-impact publications
