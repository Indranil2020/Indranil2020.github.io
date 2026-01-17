# PyMM

## Official Resources
- Homepage: https://github.com/ChenGiuseppe/PyMM
- Source Repository: https://github.com/ChenGiuseppe/PyMM
- License: MIT License

## Overview
PyMM is a Python-based software package designed for Quantum Mechanics/Molecular Mechanics (QM/MM) simulations. It specifically implements the Perturbed Matrix Method (PMM) and allows for the calculation of excited states and electronic properties of molecular systems in complex environments (such as proteins or solvents) by coupling Quantum Chemical calculations with Classical Molecular Dynamics trajectories.

**Scientific domain**: QM/MM, excited states in complex environments, perturbed matrix method
**Target user community**: Biophysicists, photochemists studying solvent effects and proteins

## Theoretical Methods
- Quantum Mechanics/Molecular Mechanics (QM/MM)
- Perturbed Matrix Method (PMM)
- Time-Dependent DFT (TDDFT) coupling
- Molecular Dynamics (MD) trajectory analysis
- Electric field effects
- Perturbed Hamiltonian diagonalization

## Capabilities (CRITICAL)
- Excited state energies in environment
- Transition properties in environment
- Trajectory post-processing
- Coupling with various MD codes
- Calculation of perturbed electronic spectra
- Conformational sampling effects

**Sources**: GitHub repository

## Key Strengths

### Perturbed Matrix Method:
- Efficient inclusion of environmental effects
- No need to re-run QC for every frame (in some approximations)
- Rigorous statistical sampling

### Python Flexibility:
- Easy to modify and extend
- Integrates with MD analysis tools (MDAnalysis)
- Flexible input handling

### Excited State Focus:
- Specifically targets properties like absorption/emission
- Perturbed transition dipoles

## Inputs & Outputs
- **Input formats**:
  - MD trajectory (GROMACS, NAMD, etc.)
  - Quantum chemical reference data (e.g. from ORCA, Gaussian)
  - Configuration files
  
- **Output data types**:
  - Perturbed energies
  - Perturbed detailed properties
  - Spectra

## Interfaces & Ecosystem
- **MD Codes**: GROMACS, NAMD (via trajectory reading)
- **QC Codes**: ORCA, Gaussian (for unperturbed reference)
- **Language**: Python

## Advanced Features

### Sampling:
- Statistical convergence of properties
- Sampling over thousands of frames
- Distribution of excitation energies

## Performance Characteristics
- **Speed**: High (post-processing)
- **Bottleneck**: Initial QC calculation
- **Scaling**: Linear with number of frames

## Computational Cost
- **Efficiency**: Very high compared to ONIOM-like QM/MM re-optimization
- **Cost**: Depends on MD length

## Limitations & Known Constraints
- **Methodology**: Restricted to PMM approximation
- **Reference**: Needs good gas-phase/reference QC calculation
- **Force Field**: Quality of MM environment matters

## Comparison with Other Codes
- **vs Chemshell**: PyMM is lighter, focused on PMM analysis
- **vs Q-Chem/Gaussian QM/MM**: PyMM is a post-processing tool for trajectories
- **Unique strength**: Perturbed Matrix Method implementation in Python

## Application Areas
- **Biochromophores**: GFP, rhodopsins
- **Solvatochromism**: Dye shifts in solvents
- **Fluctuations**: Spectral broadening due to thermal motion

## Best Practices
- **Sampling**: Ensure MD trajectory is equilibrated
- **Reference**: High-level QC for the unperturbed reference state
- **QC/MM boundary**: Define carefully
- **Convergence**: Check property distributions

## Community and Support
- Open-source MIT
- GitHub issues
- Documentation in repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ChenGiuseppe/PyMM

**Confidence**: VERIFIED - Active GitHub project

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (MIT)
- Specialized strength: QM/MM via Perturbed Matrix Method
