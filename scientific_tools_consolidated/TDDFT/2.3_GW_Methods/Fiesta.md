# Fiesta (First-principles Excited-STAte calculations)

## Official Resources
- Homepage: https://gitlab.com/marcodalessandro76/Fiesta
- Documentation: https://gitlab.com/marcodalessandro76/Fiesta/-/wikis/home
- Source Repository: https://gitlab.com/marcodalessandro76/Fiesta
- License: GNU GPL v3

## Overview
Fiesta is an open-source code for calculating electronic excitations using many-body perturbation theory (GW approximation and Bethe-Salpeter equation) starting from plane-wave DFT calculations. Developed by Marco D'Alessandro and collaborators, Fiesta focuses on efficient GW/BSE implementations with emphasis on optical properties, core-level spectroscopy, and exciton physics. It interfaces with major DFT codes (Quantum ESPRESSO, VASP) and provides comprehensive tools for excited-state calculations.

**Scientific domain**: GW approximation, BSE, optical properties, MBPT  
**Target user community**: Spectroscopy researchers, excited-state physicists, DFT users

## Theoretical Methods
- GW approximation (G₀W₀, evGW)
- Bethe-Salpeter Equation (BSE)
- Random Phase Approximation (RPA)
- Plane-wave basis
- Pseudopotentials
- Core-level excitations
- Exciton physics
- Optical properties

## Capabilities (CRITICAL)
- GW quasiparticle energies
- Accurate band gaps
- BSE optical spectra
- Exciton binding energies
- Core-level spectroscopy
- X-ray absorption (XAS)
- Optical absorption
- Interfaces with QE/VASP
- Open-source implementation
- Efficient algorithms
- Production calculations

**Sources**: Fiesta GitLab repository

## Key Strengths

### DFT Code Interfaces:
- Quantum ESPRESSO interface
- VASP interface
- Standard DFT output
- Flexible input
- Wide compatibility

### Core-Level Spectroscopy:
- XAS calculations
- Core excitations
- Element-specific
- Experimental comparison
- Comprehensive treatment

### Optical Properties:
- BSE implementation
- Exciton calculations
- Absorption spectra
- Oscillator strengths
- Accurate predictions

### Open Source:
- GNU GPL v3
- Free software
- Transparent code
- Community development
- Educational value

### Efficiency:
- Optimized algorithms
- Parallel implementation
- Production performance
- Reasonable cost
- Research quality

## Inputs & Outputs
- **Input formats**:
  - QE wavefunctions
  - VASP output
  - Fiesta input files
  - DFT ground state
  
- **Output data types**:
  - Quasiparticle energies
  - Optical spectra
  - XAS spectra
  - Exciton eigenstates
  - Band gaps

## Interfaces & Ecosystem
- **Quantum ESPRESSO**:
  - Primary interface
  - QE wavefunction input
  - Tested workflow
  - Standard integration
  
- **VASP**:
  - VASP interface
  - Alternative DFT input
  - Compatibility
  
- **Visualization**:
  - Spectral plotting
  - Analysis tools
  - Standard formats

## Workflow and Usage

### Typical Workflow:
1. Run DFT calculation (QE or VASP)
2. Prepare Fiesta input
3. Run GW calculation
4. Analyze quasiparticle energies
5. Run BSE for optical properties
6. Extract and visualize spectra

### GW Calculation:
```bash
fiesta -i gw_input
# Computes GW corrections
```

### BSE for Optics:
```bash
fiesta -i bse_input
# Solves BSE for excitons
```

## Advanced Features

### GW Implementation:
- G₀W₀ calculations
- evGW (eigenvalue SC)
- Efficient algorithms
- Plasmon-pole models
- Production quality

### BSE Capabilities:
- Electron-hole interaction
- Exciton eigenstates
- Binding energies
- Optical absorption
- Singlet excitations

### Core-Level Excitations:
- XAS implementation
- Core-hole treatment
- Element-specific spectra
- Edge calculations
- Experimental comparison

### Parallelization:
- MPI parallelization
- Efficient scaling
- Production performance
- Large systems feasible

## Performance Characteristics
- **Speed**: Good (optimized algorithms)
- **Accuracy**: Excellent for optical/XAS
- **System size**: Moderate to large
- **Scaling**: Parallel implementation
- **Typical**: Research calculations

## Computational Cost
- **GW**: Standard GW cost
- **BSE**: Moderate additional expense
- **Parallelization**: Good scaling
- **Production**: Feasible
- **Efficiency**: Competitive

## Limitations & Known Constraints
- **DFT dependency**: Requires QE or VASP
- **Learning curve**: MBPT expertise needed
- **Documentation**: Community-level
- **Support**: Developer and community
- **Platform**: Linux systems

## Comparison with Other Codes
- **vs BerkeleyGW**: Fiesta more specialized features
- **vs Yambo**: Both comprehensive GW/BSE
- **vs WEST**: Fiesta includes XAS
- **Unique strength**: Open-source, XAS capabilities, QE/VASP interfaces, core-level spectroscopy

## Application Areas

### Optical Spectroscopy:
- Absorption spectra
- Exciton physics
- Optical properties
- Experimental comparison
- Materials characterization

### X-ray Spectroscopy:
- XAS calculations
- Core-level excitations
- Element-specific
- NEXAFS
- Synchrotron experiments

### Band Gaps:
- Accurate fundamental gaps
- Quasiparticle corrections
- Semiconductor properties
- Band structure refinement

### Materials Science:
- Electronic excitations
- Optical properties
- Spectroscopy interpretation
- Materials design

## Best Practices

### DFT Preparation:
- Converged QE/VASP calculation
- Sufficient k-points
- Empty states included
- Quality ground state

### GW Convergence:
- k-point convergence
- Cutoff parameters
- Empty bands
- Frequency grid

### BSE Calculations:
- Appropriate transitions
- k-point sampling
- Exciton convergence
- Numerical parameters

### XAS Calculations:
- Core-hole treatment
- Broadening parameters
- Edge selection
- Experimental comparison

## Community and Support
- Open-source (GPL v3)
- GitLab repository
- Wiki documentation
- Developer support
- User community
- Active development

## Educational Resources
- GitLab wiki
- Example calculations
- Tutorial files
- Published papers
- User contributions

## Development
- Marco D'Alessandro (lead developer)
- Open development on GitLab
- Community contributions
- Regular updates
- Research-driven

## Research Applications
- Optical properties
- X-ray spectroscopy
- Exciton physics
- Band gap predictions
- Spectroscopy theory

## Technical Innovation

### Core-Level Focus:
- XAS implementation
- Core excitations
- Element-specific
- Comprehensive treatment

### Multiple DFT Backends:
- QE interface
- VASP interface
- Flexible framework
- Wide applicability

## Open-Source MBPT
- Free GW/BSE code
- Transparent algorithms
- Community-driven
- Educational value
- Research accessible

## Verification & Sources
**Primary sources**:
1. GitLab: https://gitlab.com/marcodalessandro76/Fiesta
2. Wiki documentation
3. M. D'Alessandro et al., publications
4. User manual

**Secondary sources**:
1. GW/BSE literature
2. XAS spectroscopy papers
3. Application studies
4. Method development

**Confidence**: VERIFIED - Active open-source code

**Verification status**: ✅ VERIFIED
- GitLab: ACCESSIBLE
- Documentation: Wiki available
- Source code: Open (GPL v3)
- Community support: Active
- Development: Regular commits
- Specialized strength: Open-source GW/BSE, core-level XAS spectroscopy, QE/VASP interfaces, optical properties, exciton calculations, production quality
