# EDMFTF

## Official Resources
- Homepage: https://hauleweb.rutgers.edu/EDMFTF/
- Documentation: https://hauleweb.rutgers.edu/EDMFTF/index.html
- Source Repository: Available via homepage registration
- License: Free for academic use (registration required)

## Overview
EDMFTF (Embedded Dynamical Mean Field Theory Functional) is a DFT+DMFT implementation developed by Kristjan Haule at Rutgers University. It provides a sophisticated interface for performing charge self-consistent DFT+DMFT calculations with advanced impurity solvers, focusing on strongly correlated materials with realistic crystal structures.

**Scientific domain**: Strongly correlated electron systems, transition metal oxides, heavy fermions, actinides  
**Target user community**: Researchers studying strongly correlated materials requiring DFT+DMFT methodology

## Theoretical Methods
- Density Functional Theory + Dynamical Mean Field Theory (DFT+DMFT)
- Charge self-consistent DFT+DMFT
- Continuous-time quantum Monte Carlo (CTQMC) impurity solver
- Hybridization expansion (CT-HYB)
- Exact diagonalization (ED) impurity solver
- Non-crossing approximation (NCA)
- One-crossing approximation (OCA)
- Hubbard I approximation
- Full Coulomb vertex implementation
- Spin-orbit coupling
- LDA+DMFT, GGA+DMFT
- Non-collinear magnetism

## Capabilities (CRITICAL)
- Charge self-consistent DFT+DMFT calculations
- Electronic structure of strongly correlated materials
- Spectral functions and DOS including correlation effects
- Magnetic properties (moments, ordering)
- Metal-insulator transitions
- Orbital ordering and occupations
- Crystal field effects in correlated systems
- Temperature-dependent properties
- Pressure-dependent studies
- Heavy fermion systems
- Actinide materials (5f electrons)
- Transition metal oxides (3d electrons)
- Rare earth systems (4f electrons)
- Momentum-resolved spectral functions
- Optical conductivity
- Thermodynamics of correlated systems

**Sources**: Official EDMFTF website (https://hauleweb.rutgers.edu/EDMFTF/), cited in 6/7 source lists

## Key Features

### Charge Self-Consistency:
- Full charge self-consistency between DFT and DMFT
- Iterative solution of DFT and DMFT equations
- Convergence acceleration techniques
- Proper treatment of charge redistribution

### Advanced Impurity Solvers:
- Multiple solver options for different regimes
- CTQMC for general multi-orbital problems
- ED for small clusters or strong coupling
- NCA/OCA for heavy fermion physics
- Solver selection based on physics

### Realistic Materials:
- Interfaces with LAPW codes (WIEN2k)
- Full crystal structure handling
- Arbitrary number of correlated orbitals
- Multiple correlated atoms per unit cell
- Spin-orbit coupling included

### Sophisticated Physics:
- Full Coulomb vertex (beyond density-density)
- Crystal field splitting
- Ligand field effects
- Covalency and hybridization
- Multi-orbital correlations

## Inputs & Outputs
- **Input formats**:
  - DMFT input file (PARAMS)
  - DFT output from WIEN2k (case.struct, etc.)
  - Wannier projections or projectors
  - Coulomb interaction parameters (U, J)
  - CTQMC parameters (beta, number of steps)
  
- **Output data types**:
  - Self-energy (Matsubara frequencies)
  - Green's functions (local and momentum-resolved)
  - Spectral functions and DOS
  - Occupation matrices
  - Magnetic moments
  - Total energy
  - Convergence data
  - Observable files for analysis

## Interfaces & Ecosystem
- **DFT code integration**:
  - Primary interface: WIEN2k (LAPW method)
  - Tight integration with WIEN2k workflow
  - Reads WIEN2k outputs directly
  
- **Impurity solver interfaces**:
  - CTQMC solver included
  - Interface to exact diagonalization
  - NCA/OCA solvers available
  
- **Wannier functions**:
  - Internal projector generation
  - Can use pre-computed Wannier functions
  - Flexible orbital selection
  
- **Analysis tools**:
  - Python scripts for post-processing
  - Spectral function analysis
  - Analytical continuation tools

## Workflow and Usage

### Typical DFT+DMFT Workflow:

1. **DFT Calculation**:
   - Run standard WIEN2k DFT calculation
   - Converge electronic structure
   - Generate necessary files

2. **Setup DMFT**:
   - Define correlated orbitals
   - Set interaction parameters (U, J)
   - Choose impurity solver
   - Set temperature and CTQMC parameters

3. **DMFT Iterations**:
   - Run DMFT self-consistency loop
   - Solve impurity problem at each iteration
   - Update DFT charge density
   - Monitor convergence

4. **Analysis**:
   - Extract spectral functions
   - Compute physical observables
   - Perform analytical continuation
   - Compare with experiments

### Parameter Selection:
- **U and J**: From constrained RPA or literature
- **Temperature**: Physical temperature or convergence parameter
- **Double counting**: Various schemes (FLL, AMF, etc.)
- **Projectors**: Atomic-like or Wannier-like

## Advanced Capabilities

### Metal-Insulator Transitions:
- Pressure or doping-induced transitions
- Temperature-driven transitions (Mott physics)
- Orbital-selective Mott transitions
- Volume collapse transitions

### Magnetic Properties:
- Paramagnetic, ferromagnetic, antiferromagnetic states
- Magnetic phase diagrams
- Spin and orbital moments
- Magnetic exchange interactions

### Spectroscopy:
- Photoemission spectroscopy (PES/ARPES) simulation
- X-ray absorption spectroscopy (XAS)
- Optical conductivity
- Momentum-resolved spectra

### Thermodynamics:
- Entropy calculations
- Specific heat
- Free energy
- Phase stability

## Computational Efficiency
- **CTQMC solver**: Most expensive component
- **Typical iteration**: Hours to days depending on system
- **Convergence**: 10-50 DMFT iterations typically
- **Parallelization**: CTQMC parallelized via MPI
- **Total calculation**: Days to weeks for production runs

## Limitations & Known Constraints
- **Computational cost**: Very expensive; CTQMC-DMFT intensive
- **WIEN2k dependency**: Requires WIEN2k for DFT part
- **Registration required**: Free but needs registration
- **Learning curve**: Very steep; requires deep understanding of DMFT
- **Statistical noise**: CTQMC introduces Monte Carlo errors
- **Analytical continuation**: MaxEnt adds systematic uncertainty
- **Parameter dependence**: Results depend on U, J, double counting
- **System size**: Limited to relatively small unit cells
- **Temperature**: Low temperatures computationally demanding
- **Documentation**: Good but assumes DMFT expertise
- **Platform**: Linux/Unix; requires proper build environment

## Application Areas

### Transition Metal Oxides:
- Cuprates (high-Tc superconductors)
- Manganites (colossal magnetoresistance)
- Vanadates, titanates
- Ruthenates

### Heavy Fermion Systems:
- Cerium and Ytterbium compounds
- Kondo lattice physics
- Valence fluctuations
- Quantum criticality

### Actinides:
- Plutonium and other 5f systems
- Volume collapse transitions
- Magnetic properties
- Electronic structure

### Iron-Based Superconductors:
- Iron pnictides and chalcogenides
- Orbital selectivity
- Magnetic order

## Comparison with Other DMFT Codes
- **vs TRIQS**: EDMFTF more focused on charge self-consistency
- **vs w2dynamics**: EDMFTF integrates with WIEN2k (LAPW)
- **vs DMFTwDFT**: Similar functionality, different implementations
- **Unique strength**: Sophisticated charge self-consistency, Haule's expertise

## Best Practices

### Convergence:
- Start with Hubbard I or non-self-consistent
- Gradually increase CTQMC accuracy
- Monitor charge self-consistency carefully
- Use charge mixing for stability

### Parameter Choice:
- Validate U, J from constrained calculations
- Test double counting scheme sensitivity
- Consider multiple temperature points
- Check projector dependence

### Verification:
- Compare with experiments (PES, XAS, etc.)
- Check sum rules and conservation laws
- Validate against known limits
- Perform convergence studies

## Verification & Sources
**Primary sources**:
1. Official website: https://hauleweb.rutgers.edu/EDMFTF/
2. Documentation: https://hauleweb.rutgers.edu/EDMFTF/index.html
3. K. Haule et al., Phys. Rev. B 81, 195107 (2010) - Charge self-consistent DFT+DMFT
4. K. Haule, Phys. Rev. B 75, 155113 (2007) - Quantum Monte Carlo impurity solver
5. K. Haule and G. Kotliar, New J. Phys. 11, 025021 (2009) - Coherence-incoherence crossover
6. P. Werner and A. J. Millis, Phys. Rev. B 74, 155107 (2006) - CT-HYB algorithm

**Secondary sources**:
1. EDMFTF tutorials and workshops
2. Published DFT+DMFT studies from Haule group
3. Rutgers strongly correlated materials research
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Software: Free for academics (registration required)
- Community support: Active (Haule group, email support)
- Academic citations: >300 (method and application papers)
- Active use: Standard for charge self-consistent DFT+DMFT
- Benchmark validation: Extensive comparisons with experiments
- Developer: K. Haule (leading DMFT expert at Rutgers)
