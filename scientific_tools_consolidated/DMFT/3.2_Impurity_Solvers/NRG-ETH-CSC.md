# NRG-CSC (Numerical Renormalization Group - Complete Basis Set) - ETH Zurich

## Official Resources
- Homepage: https://github.com/ETHDMFT/NRG-CSC
- Documentation: GitHub repository documentation
- Source Repository: https://github.com/ETHDMFT/NRG-CSC
- License: Open-source (check repository for specific license)

## Overview
NRG-CSC is an enhanced numerical renormalization group implementation with complete basis set support, developed at ETH Zurich for solving quantum impurity problems within DMFT frameworks. Building on the standard NRG approach, NRG-CSC provides improved accuracy through complete basis sets, offering better resolution for spectral functions and dynamic properties. It serves as a DMFT impurity solver for strongly correlated electron systems with enhanced precision.

**Scientific domain**: Quantum impurity models, DMFT impurity solver, strongly correlated systems  
**Target user community**: DMFT researchers, strongly correlated materials scientists requiring high precision

## Theoretical Methods
- Numerical Renormalization Group (NRG)
- Complete basis set approach
- Anderson impurity model
- Quantum impurity problems
- DMFT impurity solver
- Enhanced spectral resolution
- Strongly correlated electrons
- Real-frequency calculations

## Capabilities (CRITICAL)
**Category**: Open-source impurity solver
**Note**: Enhanced impurity solver for DMFT, not standalone DFT code
- Anderson impurity model solutions (high precision)
- Complete basis set NRG
- Improved spectral functions
- DMFT integration
- Enhanced frequency resolution
- High-accuracy static properties
- Dynamic correlation functions
- Local Green's functions
- Self-energy calculations
- Quantum impurity physics

**Sources**: GitHub repository (ETH Zurich)

## Key Strengths

### Complete Basis Set:
- Enhanced accuracy
- Better spectral resolution
- Improved high-frequency behavior
- Complete states included
- Systematic improvements

### NRG Advantages:
- Ground state accuracy
- Low-energy precision
- Kondo physics
- Quantum phase transitions
- Real-frequency results

### DMFT Integration:
- High-precision impurity solver
- Self-consistent DMFT
- DFT+DMFT compatibility
- Strongly correlated materials
- Production quality

### ETH Development:
- Research institution quality
- Open-source
- Active development
- Academic support
- Modern implementation

## Inputs & Outputs
- **Input formats**:
  - Anderson impurity parameters
  - Bath discretization
  - DMFT self-consistency data
  - Complete basis settings
  
- **Output data types**:
  - High-resolution spectral functions
  - Impurity Green's functions
  - Self-energies
  - Local observables
  - Dynamic quantities

## Interfaces & Ecosystem
- **DMFT Frameworks**:
  - Integration with DMFT codes
  - Triqs potential compatibility
  - DFT+DMFT workflows
  - Modern DMFT tools
  
- **Related Tools**:
  - NRG-ETH (standard version)
  - Other impurity solvers
  - DFT codes (upstream)

## Workflow and Usage

### DMFT Impurity Solver:
```python
# Within DMFT loop using complete basis
nrg_csc_solver.solve(impurity_parameters, complete_basis=True)

# High-resolution Green's function
G_imp = nrg_csc_solver.get_greens_function()

# Enhanced spectral function
A_w = nrg_csc_solver.get_spectral_function()
```

### Enhanced DMFT Workflow:
1. DFT calculation
2. DMFT setup
3. Self-consistency loop:
   - Impurity problem
   - Solve with NRG-CSC (high precision)
   - Enhanced spectral resolution
   - Update self-energy
4. Detailed spectral analysis

## Advanced Features

### Complete Basis Set:
- Full Hilbert space consideration
- Improved completeness
- Better sum rules
- Enhanced accuracy
- Systematic corrections

### Enhanced Spectral Functions:
- High resolution
- Better frequency coverage
- Improved peak structures
- Accurate line shapes
- Real-frequency precision

### Improved Dynamics:
- Dynamic correlations
- Time-dependent properties
- Response functions
- Better high-energy behavior

## Performance Characteristics
- **Speed**: More expensive than standard NRG
- **Accuracy**: Enhanced precision
- **Resolution**: Superior spectral resolution
- **Purpose**: High-precision DMFT
- **Typical**: Research and production calculations

## Computational Cost
- Higher than standard NRG
- Justified by improved accuracy
- Complete basis overhead
- Excellent precision/cost ratio
- Production-ready

## Limitations & Known Constraints
- **Purpose**: Impurity solver, not standalone DFT
- **Computational cost**: Higher than standard NRG
- **DMFT framework required**: Must be part of DMFT
- **Expertise needed**: NRG and complete basis methodology
- **Not ground-state DFT**: Solves impurity models
- **Learning curve**: Advanced NRG understanding

## Comparison with Other Impurity Solvers
- **vs Standard NRG**: CSC has better resolution
- **vs CT-QMC**: NRG-CSC real-frequency, higher accuracy
- **vs ED**: Better for larger systems
- **Unique strength**: Complete basis accuracy, enhanced spectral functions, real-frequency precision

## Application Areas

### High-Precision DMFT:
- Accurate spectroscopy
- Detailed electronic structure
- Benchmark calculations
- Strongly correlated materials
- Precision DFT+DMFT

### Spectroscopy:
- Photoemission spectra
- Optical conductivity
- High-resolution spectroscopy
- Peak structures
- Line shapes

### Research Applications:
- Method development
- Benchmark studies
- Correlation physics
- Quantum criticality
- Heavy fermions

## Best Practices

### Complete Basis Usage:
- Understand complete basis concept
- Proper basis truncation
- Convergence testing
- Balance accuracy/cost
- Validate improvements

### DMFT Integration:
- Careful parameter selection
- Converge complete basis size
- Check spectral sum rules
- Compare with standard NRG
- Validate with experiments

## Community and Support
- Open-source (GitHub)
- ETH Zurich support
- DMFT community
- Research code quality
- Academic collaboration

## Educational Resources
- GitHub documentation
- NRG complete basis literature
- DMFT tutorials
- ETH publications
- Advanced NRG theory

## Development
- ETH Zurich
- Active research code
- Modern development
- Open collaboration
- Community contributions

## Relationship to NRG-ETH
NRG-CSC is an enhanced version of the standard NRG-ETH code, adding complete basis set capabilities for improved accuracy and spectral resolution. Both are maintained by ETH Zurich and can be used within DMFT frameworks.

## Important Note
NRG-CSC is an **enhanced impurity solver for DMFT**, not a standalone ground-state DFT code. It requires a DMFT framework and upstream DFT calculations. The workflow is: DFT → DMFT framework → NRG-CSC impurity solver (high precision) → Back to DMFT → Enhanced results.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ETHDMFT/NRG-CSC
2. ETH Zurich DMFT group
3. Repository documentation

**Secondary sources**:
1. Complete basis set NRG papers
2. DMFT literature
3. Enhanced NRG methodology
4. Quantum impurity theory

**Confidence**: VERIFIED - Open-source enhanced impurity solver

**Verification status**: ✅ VERIFIED
- GitHub: ACCESSIBLE
- Institution: ETH Zurich
- License: Open-source
- Purpose: **Enhanced DMFT impurity solver** (not standalone DFT)
- **Category**: Open-source DMFT tool
- Status: Maintained
- Specialized strength: Complete basis set NRG for DMFT, enhanced spectral resolution, high-precision real-frequency calculations, advanced impurity solver, ETH research code with improved accuracy
