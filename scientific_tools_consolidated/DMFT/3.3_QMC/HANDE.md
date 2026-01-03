# HANDE (Highly Accurate N-DEterminant)

## Official Resources
- Homepage: https://hande.readthedocs.io/
- Documentation: https://hande.readthedocs.io/
- Source Repository: https://github.com/hande-qmc/hande
- License: LGPL v2.1

## Overview
HANDE is a modern, efficient implementation of FCIQMC (Full Configuration Interaction Quantum Monte Carlo) and related stochastic quantum chemistry methods. Developed as a community code with contributions from multiple institutions, HANDE provides production-quality implementations of FCIQMC, CCMC (Coupled Cluster Monte Carlo), and DMQMC (Density Matrix QMC) with emphasis on code quality, documentation, and ease of use. The code is designed for both research applications and method development.

**Scientific domain**: FCIQMC, stochastic quantum chemistry, many-body methods  
**Target user community**: Quantum chemists, method developers, FCIQMC users

## Theoretical Methods
- Full Configuration Interaction QMC (FCIQMC)
- Coupled Cluster Monte Carlo (CCMC)
- Density Matrix QMC (DMQMC)
- Initiator approximation
- Semi-stochastic adaptations
- Finite-temperature methods
- Excited states
- Exact diagonalization

## Capabilities (CRITICAL)
**Category**: Open-source FCIQMC/stochastic quantum chemistry
- FCIQMC implementation
- CCMC variants
- DMQMC (finite-T)
- Semi-stochastic methods
- Initiator FCIQMC
- Molecules and lattices
- Ground and excited states
- Real and complex integrals
- Symmetry exploitation
- MPI + OpenMP parallelization
- HDF5 output
- Python interface (pyHANDE)
- Production quality

**Sources**: Official documentation, GitHub, publications

## Key Strengths

### Modern Implementation:
- Clean Fortran code
- Well-documented
- Unit tested
- Community-driven
- Best practices

### Comprehensive Methods:
- FCIQMC, CCMC, DMQMC
- Multiple approaches
- Flexibility
- Research capabilities
- Production quality

### Python Integration:
- pyHANDE interface
- Workflow automation
- Analysis tools
- Visualization
- User-friendly

### Code Quality:
- Extensive testing
- Clear documentation
- Active development
- Issue tracking
- Community support

## Inputs & Outputs
- **Input formats**:
  - Lua input scripts
  - FCIDUMP integrals
  - Model Hamiltonians
  - Configuration files
  
- **Output data types**:
  - Energies and observables
  - HDF5 archives
  - Sampling data
  - Analysis-ready formats
  - pyHANDE objects

## Interfaces & Ecosystem

### Quantum Chemistry:
- FCIDUMP format
- Molpro
- PySCF
- GAMESS
- Dalton

### Python Tools:
- pyHANDE analysis
- Jupyter workflows
- Visualization
- Data management

## Workflow and Usage

### Installation:
```bash
# Clone repository
git clone https://github.com/hande-qmc/hande.git
cd hande
mkdir build && cd build
cmake ..
make -j8
```

### Lua Input (fciqmc.lua):
```lua
sys = hubbard_k {
    electrons = 8,
    lattice = { {4} },
    ms = 0,
    U = 1.3,
    t = 1.0,
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.01,
        rng_seed = 7,
        init_pop = 10,
        mc_cycles = 10,
        nreports = 100,
        target_population = 1e6,
        state_size = -500,
        spawned_state_size = -100,
    },
}
```

### Run HANDE:
```bash
# MPI + OpenMP
export OMP_NUM_THREADS=4
mpirun -n 4 hande.x fciqmc.lua > output.out
```

### Python Analysis:
```python
import pyhande

# Load HANDE output
data = pyhande.extract.extract_data('output.out')

# Analyze results
(data_len, reblock_data, covariance) = pyhande.analysis.reblock(data[0])

# Plot
pyhande.analysis.plot_reblocking(reblock_data)
```

## Advanced Features

### Semi-Stochastic:
- Deterministic space
- Stochastic remainder
- Improved efficiency
- Reduced noise
- Larger systems

### Symmetry:
- Point group symmetries
- Momentum conservation
- Spin symmetries
- Computational efficiency

### Finite Temperature:
- DMQMC method
- Thermal properties
- Temperature-dependent
- Phase transitions
- Statistical mechanics

### Real-Time:
- Time-dependent methods
- Dynamics (development)
- Excited states
- Response properties

## Performance Characteristics
- **Speed**: Efficient MPI+OpenMP
- **Accuracy**: Numerically exact (converged)
- **Scalability**: Good parallel scaling
- **System size**: Moderate (CI space)
- **Purpose**: Strongly correlated, benchmarks

## Computational Cost
- Walker population dependent
- CI space scaling
- Expensive but exact
- HPC suitable
- Production capable

## Limitations & Known Constraints
- **CI space**: Basis set limited
- **Computational cost**: Expensive
- **Sign problem**: Fermions
- **System size**: Moderate
- **HPC recommended**: Production calculations

## Comparison with Other FCIQMC Codes
- **vs NECI**: HANDE similar methods, modern codebase
- **vs Traditional FCI**: HANDE stochastic, larger systems
- **Unique strength**: Modern implementation, Python integration, community code, documentation, code quality

## Application Areas

### Strongly Correlated Chemistry:
- Molecules
- Transition metals
- Multi-reference
- Bond breaking
- Strong correlation

### Benchmarks:
- Exact energies
- Method validation
- Reference results
- Correlation energies
- Chemical accuracy

### Lattice Models:
- Hubbard model
- t-J model
- Quantum magnetism
- Model systems

### Method Development:
- FCIQMC research
- Algorithm innovation
- Stochastic methods
- Semi-stochastic

## Best Practices

### Input Scripts:
- Lua flexibility
- Clear syntax
- Modular approach
- Documentation

### Analysis:
- Use pyHANDE
- Reblocking analysis
- Error estimation
- Statistical rigor

### Production:
- Sufficient walkers
- Convergence testing
- Multiple runs
- HPC resources

## Community and Support
- Open-source (LGPL v2.1)
- Multi-institutional
- GitHub repository
- Active development
- Comprehensive docs
- Community-driven
- Issue tracking

## Educational Resources
- Excellent documentation
- Tutorials
- Example inputs
- pyHANDE guides
- FCIQMC literature
- API reference

## Development
- Community collaboration
- Multiple institutions
- Active development
- Modern practices
- Regular releases
- User feedback

## Research Impact
HANDE provides accessible, well-documented FCIQMC implementations, enabling researchers to apply cutting-edge stochastic quantum chemistry methods with confidence in code quality and results.

## Verification & Sources
**Primary sources**:
1. Homepage: https://hande.readthedocs.io/
2. GitHub: https://github.com/hande-qmc/hande
3. Publications: J. Chem. Theory Comput. 15, 3 (2019)

**Secondary sources**:
1. FCIQMC literature
2. User publications
3. Quantum chemistry papers

**Confidence**: VERIFIED - Modern FCIQMC code

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: LGPL v2.1 (open-source)
- **Category**: Open-source FCIQMC code
- Status: Actively developed
- Community: Multi-institutional
- Specialized strength: Modern FCIQMC/CCMC/DMQMC implementation, excellent documentation, Python integration (pyHANDE), code quality emphasis, community-driven, production quality, semi-stochastic methods, Lua input, HDF5 output, comprehensive testing, user-friendly
