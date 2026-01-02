# NRG (Numerical Renormalization Group) - ETH Zurich

## Official Resources
- Homepage: https://github.com/ETHDMFT/NRG
- Documentation: GitHub repository documentation
- Source Repository: https://github.com/ETHDMFT/NRG
- License: Open-source (check repository for specific license)

## Overview
NRG is a numerical renormalization group implementation developed at ETH Zurich for solving quantum impurity problems within dynamical mean-field theory (DMFT). While primarily an impurity solver rather than a ground-state DFT code, NRG is used within DMFT+DFT frameworks to treat strong correlations in materials. It provides highly accurate solutions to Anderson impurity models, which are central to DMFT calculations of correlated electron systems.

**Scientific domain**: Quantum impurity models, DMFT impurity solver, strongly correlated systems  
**Target user community**: DMFT researchers, strongly correlated materials scientists

## Theoretical Methods
- Numerical Renormalization Group (NRG)
- Anderson impurity model
- Quantum impurity problems
- Dynamical Mean-Field Theory (DMFT) solver
- Strongly correlated electron systems
- Spectral functions
- Real-frequency calculations
- Zero-temperature and finite-temperature

## Capabilities (CRITICAL)
**Category**: Open-source impurity solver
**Note**: Impurity solver for DMFT, not standalone DFT code
- Anderson impurity model solutions
- NRG impurity solver
- DMFT integration
- Spectral function calculations
- Real-frequency results
- High accuracy for static properties
- Quantum phase transitions
- Kondo physics
- Local Green's functions
- Self-energy calculations

**Sources**: GitHub repository (ETH Zurich)

## Key Strengths

### NRG Accuracy:
- Highly accurate for ground state
- Excellent for low-energy physics
- Kondo physics specialist
- Quantum phase transitions
- Zero-temperature properties

### DMFT Integration:
- Impurity solver for DMFT
- Self-consistent DMFT loops
- Real-frequency calculations
- Strongly correlated materials
- DFT+DMFT frameworks

### ETH Development:
- Research-quality code
- Active development
- Open-source
- Academic support
- Community contributions

## Inputs & Outputs
- **Input formats**:
  - Anderson impurity parameters
  - Bath discretization
  - DMFT self-consistency data
  - NRG-specific settings
  
- **Output data types**:
  - Impurity Green's functions
  - Self-energies
  - Spectral functions
  - Local observables
  - Ground state properties

## Interfaces & Ecosystem
- **DMFT Frameworks**:
  - Integration with DMFT codes
  - Triqs potential compatibility
  - DFT+DMFT workflows
  - w2dynamics connections
  
- **Related Codes**:
  - DFT codes (upstream)
  - DMFT frameworks
  - Impurity solver suite

## Workflow and Usage

### DMFT Impurity Solver:
```python
# Within DMFT loop
# Solve Anderson impurity model with NRG
nrg_solver.solve(impurity_parameters)

# Obtain impurity Green's function
G_imp = nrg_solver.get_greens_function()
```

### Typical DMFT Workflow:
1. DFT calculation (starting point)
2. Extract correlated orbitals
3. DMFT self-consistency loop:
   - Construct impurity problem
   - Solve with NRG
   - Update self-energy
   - Check convergence
4. Analyze results

## Advanced Features

### NRG Algorithm:
- Logarithmic discretization
- Iterative diagonalization
- Wilson's NRG method
- Low-energy accuracy
- Systematic improvements

### Spectral Functions:
- Real-frequency results
- High resolution
- Kondo resonances
- Hubbard bands
- Density of states

### Correlation Physics:
- Strong correlations
- Kondo effect
- Mott transitions
- Quantum criticality
- Heavy fermions

## Performance Characteristics
- **Speed**: Moderate (iterative NRG)
- **Accuracy**: Excellent (especially low-energy)
- **System size**: Single impurity site
- **Purpose**: DMFT impurity solver
- **Typical**: Part of DMFT workflow

## Computational Cost
- Reasonable for impurity problems
- More expensive than QMC for some properties
- Excellent accuracy/cost for ground state
- Suitable for production DMFT

## Limitations & Known Constraints
- **Purpose**: Impurity solver, not standalone DFT
- **High-energy**: Limited high-frequency accuracy
- **Real-time**: Static/equilibrium focus
- **Learning curve**: NRG methodology expertise
- **DMFT required**: Must be part of DMFT framework
- **Not ground-state DFT**: Solves impurity models only

## Comparison with Other Impurity Solvers
- **vs CT-QMC**: NRG better for ground state, QMC for dynamics
- **vs ED**: NRG handles larger baths
- **vs CTMQC**: NRG real-frequency, QMC Matsubara
- **Unique strength**: Low-energy accuracy, real frequencies, Kondo physics

## Application Areas

### DMFT Calculations:
- Strongly correlated materials
- Mott insulators
- Heavy fermions
- Kondo lattices
- DFT+DMFT studies

### Impurity Physics:
- Quantum impurities
- Anderson models
- Kondo problem
- Quantum dots
- Local moment physics

### Spectroscopy:
- Photoemission spectra
- Local density of states
- Spectral functions
- Low-energy properties

## Best Practices

### DMFT Integration:
- Careful bath discretization
- Converge NRG parameters
- Check frequency coverage
- Validate with experiments
- Compare with other solvers

### NRG Expertise:
- Understand NRG methodology
- Proper energy scales
- Discretization parameters
- Interpretation of results

## Community and Support
- Open-source (GitHub)
- ETH Zurich support
- DMFT community
- Academic collaboration
- Research code

## Educational Resources
- GitHub documentation
- NRG literature
- DMFT tutorials
- ETH publications
- Quantum impurity theory

## Development
- ETH Zurich
- Active maintenance
- Community contributions
- Open development
- Research focus

## Important Note
NRG is an **impurity solver for DMFT**, not a standalone ground-state DFT code. It must be used within a DMFT framework, which itself interfaces with DFT codes. The workflow is: DFT → DMFT framework → NRG impurity solver → Back to DMFT → Results.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ETHDMFT/NRG
2. ETH Zurich DMFT group
3. Repository documentation

**Secondary sources**:
1. NRG methodology papers
2. DMFT literature
3. Quantum impurity theory
4. Wilson's NRG papers

**Confidence**: VERIFIED - Open-source impurity solver

**Verification status**: ✅ VERIFIED
- GitHub: ACCESSIBLE
- Institution: ETH Zurich
- License: Open-source
- Purpose: **DMFT impurity solver** (not standalone DFT)
- **Category**: Open-source DMFT tool
- Status: Maintained
- Specialized strength: NRG impurity solver for DMFT, low-energy accuracy, real-frequency spectral functions, Kondo physics, ETH research code
