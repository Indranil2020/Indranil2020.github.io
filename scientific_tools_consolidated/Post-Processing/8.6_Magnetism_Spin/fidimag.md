# fidimag

## Official Resources
- Source Repository: https://github.com/computationalmodelling/fidimag
- Documentation: https://fidimag.readthedocs.io/
- License: BSD License

## Overview
**fidimag** (Finite DIfference microMAGnetic code) is a Python/Cython/C package for finite-difference micromagnetic and atomistic simulations. It supports both continuum micromagnetic and atomistic spin models, making it suitable for multiscale magnetic simulations.

**Scientific domain**: Micromagnetic and atomistic spin simulation  
**Target user community**: Researchers needing both micromagnetic and atomistic spin simulations in a single code

## Theoretical Methods
- Landau-Lifshitz-Gilbert (LLG) equation
- Finite-difference method
- Micromagnetic continuum model
- Atomistic Heisenberg model
- Monte Carlo simulation
- Nudged elastic band (NEB) for energy barriers
- Exchange, anisotropy, Zeeman, demagnetization energies
- Dzyaloshinskii-Moriya interaction

## Capabilities (CRITICAL)
- Micromagnetic simulation (continuum)
- Atomistic spin simulation
- Monte Carlo simulation
- Energy barrier calculation (NEB)
- Domain wall dynamics
- Skyrmion simulation
- Hysteresis loops
- DMI support
- Python scripting interface

**Sources**: GitHub repository, J. Open Res. Software 6, 22 (2018)

## Key Strengths

### Dual-Scale Simulation:
- Continuum micromagnetic model
- Atomistic Heisenberg model
- Seamless switching between models
- Multiscale capability

### NEB for Energy Barriers:
- Nudged elastic band method
- Transition path calculation
- Energy barrier determination
- Switching field estimation

### Python Interface:
- Full Python scripting
- Jupyter notebook compatible
- Easy post-processing
- Extensible framework

## Inputs & Outputs
- **Input formats**:
  - Python scripts
  - Material parameters
  - Mesh specifications
  
- **Output data types**:
  - Magnetization fields
  - Energy vs time
  - NEB paths and barriers
  - VTK output for visualization

## Interfaces & Ecosystem
- **Python**: Primary interface
- **Cython/C**: Performance-critical code
- **NumPy**: Data handling
- **Matplotlib**: Visualization

## Performance Characteristics
- **Speed**: Moderate (Cython optimized)
- **Accuracy**: Good (validated)
- **System size**: Hundreds of thousands of spins
- **Parallelization**: Limited

## Computational Cost
- **Small systems**: Minutes
- **Large systems**: Hours
- **NEB**: Hours (multiple images)
- **Typical**: Moderate

## Limitations & Known Constraints
- **No GPU**: CPU-only
- **Limited parallelization**: Mostly serial
- **Finite differences only**: No FEM
- **Community**: Smaller than OOMMF
- **Documentation**: Could be more extensive

## Comparison with Other Codes
- **vs OOMMF**: fidimag has atomistic + NEB, OOMMF is NIST standard
- **vs Spirit**: fidimag is Python, Spirit is C++ with GUI
- **vs VAMPIRE**: fidimag has NEB, VAMPIRE is more established
- **Unique strength**: Dual micromagnetic+atomistic simulation, NEB energy barriers, Python interface

## Application Areas

### Domain Walls:
- Domain wall profiles
- Pinning and depinning
- Current-driven motion
- Walker breakdown

### Skyrmions:
- Skyrmion creation
- Skyrmion Hall effect
- Skyrmion stability
- DMI-driven chirality

### Energy Barriers:
- Switching barriers
- Coercivity estimation
- Thermal stability
- Transition paths

### Multiscale:
- Atomistic-to-continuum
- Local atomistic regions
- Hybrid simulations
- Parameter extraction

## Best Practices

### Model Selection:
- Use atomistic for small/nano systems
- Use continuum for larger systems
- Compare both models for validation
- Use NEB for energy barriers

### NEB Calculations:
- Use sufficient images
- Converge spring constants
- Validate endpoint structures
- Check path smoothness

## Community and Support
- Open source (BSD)
- Developed at University of Southampton
- Published in J. Open Res. Software
- ReadTheDocs documentation
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/computationalmodelling/fidimag
2. M.-A. Bisotti et al., J. Open Res. Software 6, 22 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- Published methodology: JORS
- Active development: Maintained
- Specialized strength: Dual micromagnetic+atomistic simulation, NEB energy barriers, Python interface
