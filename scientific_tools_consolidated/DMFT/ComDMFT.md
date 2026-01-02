# ComDMFT

## Official Resources
- Homepage: https://github.com/ComDMFT/ComDMFT
- Documentation: https://github.com/ComDMFT/ComDMFT/wiki
- Source Repository: https://github.com/ComDMFT/ComDMFT
- License: GNU General Public License v3.0

## Overview
ComDMFT (Combination of DMFT codes) is an integrated DFT+DMFT package that combines multiple DFT codes with sophisticated DMFT impurity solvers. Developed by the Comscope group, it provides a unified framework for charge self-consistent DFT+DMFT calculations with emphasis on transition metal systems and strongly correlated materials.

**Scientific domain**: Strongly correlated materials, DFT+DMFT methodology, transition metal compounds  
**Target user community**: Researchers performing DFT+DMFT calculations on correlated electron systems

## Theoretical Methods
- DFT+DMFT (charge self-consistent)
- LDA+DMFT, GGA+DMFT
- Continuous-time quantum Monte Carlo (CTQMC)
- CT-HYB (hybridization expansion)
- Exact diagonalization (ED)
- Hubbard I approximation
- Spin-orbit coupling treatment
- Non-collinear magnetism
- GW+DMFT extensions
- Dual fermion approaches
- Cluster DMFT extensions

## Capabilities (CRITICAL)
- Charge self-consistent DFT+DMFT
- Multiple DFT code backends (VASP, Wannier90, etc.)
- Multi-orbital strongly correlated systems
- Spectral functions with many-body effects
- Magnetic ordering calculations
- Metal-insulator transitions
- Orbital-selective correlations
- Temperature-dependent properties
- Momentum-resolved spectroscopy
- Optical properties with correlations
- Crystal field effects
- Covalency and ligand field
- Integration with ComCTQMC solver
- Flexible impurity solver selection
- Automated workflow management
- Python-based framework

**Sources**: Official ComDMFT repository (https://github.com/ComDMFT/ComDMFT), cited in 6/7 source lists

## Key Features

### Modular Architecture:
- Separation of DFT and DMFT components
- Multiple DFT code support
- Pluggable impurity solvers
- Flexible workflow customization
- Python-based scripting

### DFT Code Integration:
- VASP interface (primary)
- Wannier90 for downfolding
- Support for other DFT codes via adapters
- Automated interface generation
- Seamless data exchange

### Impurity Solver Support:
- ComCTQMC (included CT-HYB solver)
- Interface to external solvers
- Solver parameter optimization
- Multiple solver comparison

### Charge Self-Consistency:
- Full DFT-DMFT charge self-consistency
- Iterative convergence algorithms
- Mixing schemes for stability
- Convergence monitoring and diagnostics

## Inputs & Outputs
- **Input formats**:
  - Python configuration scripts
  - VASP POSCAR and output files
  - Wannier90 outputs
  - DMFT parameter files
  - Interaction parameters (U, J matrices)
  
- **Output data types**:
  - Self-energy functions
  - Green's functions (local and k-resolved)
  - Spectral functions and DOS
  - Charge densities
  - Magnetic moments
  - Occupation matrices
  - Convergence histories
  - HDF5 data files

## Interfaces & Ecosystem
- **DFT interfaces**:
  - VASP (primary backend)
  - Wannier90 for projections
  - Interface layer for other codes
  
- **Impurity solvers**:
  - ComCTQMC (built-in)
  - TRIQS solvers (compatible)
  - Custom solver interfaces
  
- **Analysis tools**:
  - Python post-processing scripts
  - Spectral function analysis
  - MaxEnt analytical continuation
  - Visualization utilities
  
- **Workflow management**:
  - Python-based workflow scripts
  - Automated job submission
  - Checkpoint and restart

## Workflow and Usage

### Typical DFT+DMFT Workflow:

1. **DFT Preparation**:
   - Run VASP DFT calculation
   - Generate Wannier functions with Wannier90
   - Define correlated orbitals

2. **DMFT Setup**:
   - Configure ComDMFT parameters
   - Set interaction parameters (U, J)
   - Choose impurity solver and settings
   - Define temperature and frequency grids

3. **Self-Consistent Calculation**:
   - Initialize DMFT loop
   - Iterate DFT and DMFT steps
   - Solve impurity problem each iteration
   - Update charge density
   - Monitor convergence

4. **Post-Processing**:
   - Extract spectral functions
   - Perform analytical continuation
   - Calculate observables
   - Compare with experiments

### Python Scripting Example:
```python
# Example workflow structure
import comdmft

# Initialize calculation
calc = comdmft.DMFTCalculation(
    dft_code='vasp',
    structure='POSCAR',
    correlated_orbitals=['d']
)

# Set parameters
calc.set_hubbard_u(U=5.0, J=0.7)
calc.set_temperature(T=300)

# Run DFT+DMFT
calc.run_self_consistent(max_iter=50)

# Analyze results
spectrum = calc.get_spectral_function()
```

## Advanced Features

### Orbital Selection:
- Flexible definition of correlated subspace
- Wannier-based downfolding
- Energy window selection
- Orbital character analysis

### Interaction Parameters:
- Full Coulomb interaction matrix
- Kanamori parametrization
- Slater integrals
- Constrained RPA calculations
- Screening effects

### Double Counting:
- Multiple schemes (FLL, AMF, etc.)
- Orbital-dependent corrections
- Self-consistent determination

### Convergence Acceleration:
- Charge mixing algorithms
- DIIS extrapolation
- Anderson mixing
- Adaptive damping

## Computational Aspects

### Performance:
- CTQMC solver: computationally intensive
- Typical iteration: hours to days
- Full calculation: days to weeks
- MPI parallelization for solver
- Python overhead minimal

### Memory Requirements:
- Moderate for DMFT framework
- Solver-dependent memory usage
- CTQMC most memory-intensive
- HDF5 for efficient data storage

### Scalability:
- Good parallel scaling for CTQMC
- Multiple k-point parallelization
- Frequency parallelization possible

## Limitations & Known Constraints
- **Computational cost**: DFT+DMFT very expensive
- **VASP focus**: Primarily designed for VASP
- **Learning curve**: Steep; requires DFT+DMFT knowledge
- **Python dependency**: Requires Python environment
- **CTQMC limitations**: Sign problem, statistical noise
- **Analytical continuation**: MaxEnt uncertainties
- **Parameter dependence**: Results sensitive to U, J, double counting
- **System size**: Limited by DFT and DMFT costs
- **Documentation**: Good but assumes expertise
- **Platform**: Linux/Unix

## Application Areas

### Transition Metal Oxides:
- Cuprates and nickelates
- Cobaltites and manganites
- Vanadates and titanates
- Correlation-driven phenomena

### Strongly Correlated Materials:
- Mott insulators
- Heavy fermion systems
- Orbital ordering
- Magnetic phase transitions

### Spectroscopy Simulation:
- Photoemission spectroscopy (ARPES)
- X-ray absorption (XAS/XMCD)
- Optical conductivity
- Transport properties

### Materials Design:
- Screening correlated materials
- Property prediction
- Phase diagram exploration

## Comparison with Other DMFT Frameworks
- **vs TRIQS**: ComDMFT more focused on VASP integration
- **vs EDMFTF**: Different DFT backend, similar capabilities
- **vs w2dynamics**: ComDMFT full DFT+DMFT framework vs solver only
- **Unique strength**: Python-based, modular architecture, VASP focus

## Best Practices

### Getting Started:
1. Start with non-self-consistent calculation
2. Validate with known systems
3. Carefully converge parameters
4. Document all settings

### Parameter Selection:
- Validate U and J from literature or cRPA
- Test multiple double counting schemes
- Explore temperature dependence
- Check frequency grid convergence

### Convergence:
- Use mixing for stability
- Monitor all observables
- Save checkpoints frequently
- Compare different solver settings

### Validation:
- Compare with experiments
- Check sum rules
- Validate limiting cases
- Perform systematic studies

## Community and Development
- Active development on GitHub
- Open-source contributions welcome
- Issue tracking and bug reports
- Community discussions
- Regular updates and improvements

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/ComDMFT/ComDMFT
2. Documentation: https://github.com/ComDMFT/ComDMFT/wiki
3. ComCTQMC solver: https://github.com/ComDMFT/ComCTQMC
4. Related publications from Comscope group

**Secondary sources**:
1. ComDMFT tutorials and examples
2. Workshop materials
3. Published DFT+DMFT studies using ComDMFT
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (Wiki)
- Source code: OPEN (GitHub, GPL v3)
- Community support: Active (GitHub issues, discussions)
- Academic citations: Growing user base
- Active development: Regular commits and updates
- Integration: Well-integrated with VASP and Wannier90
- Framework: Modern Python-based architecture
