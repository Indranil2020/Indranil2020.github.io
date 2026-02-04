# EDRIXS

## Official Resources
- Homepage: https://nsls-ii.github.io/edrixs/
- GitHub: https://github.com/NSLS-II/edrixs
- Documentation: https://nsls-ii.github.io/edrixs/
- Publication: Y.L. Wang et al., Comput. Phys. Commun. 243, 151 (2019)
- License: GNU General Public License v3.0

## Overview
EDRIXS is an open-source toolkit for simulating resonant inelastic X-ray scattering (RIXS) spectra of strongly correlated materials. It uses exact diagonalization of model Hamiltonians including crystal field, spin-orbit coupling, and Coulomb interactions to calculate XAS and RIXS spectra at transition metal and rare earth edges.

**Scientific domain**: RIXS spectroscopy, strongly correlated materials, transition metal oxides
**Target user community**: Researchers studying correlated electron systems with RIXS

## Theoretical Methods
- Exact diagonalization of many-body Hamiltonians
- Crystal field theory
- Atomic multiplet calculations
- Spin-orbit coupling
- Coulomb interaction (Slater integrals)
- RIXS cross-section calculation

## Capabilities (CRITICAL)
- **RIXS Simulation**: Full RIXS spectra calculation
- **XAS Calculation**: X-ray absorption at various edges
- **Exact Diagonalization**: Many-body eigenstates
- **Crystal Field**: Arbitrary symmetry support
- **Spin-Orbit**: Full relativistic treatment
- **Python API**: Scriptable interface
- **Parallel**: MPI parallelization

**Sources**: EDRIXS documentation, CPC publication

## Key Strengths

### Exact Methods:
- Full multiplet theory
- No approximations on correlations
- Correct symmetry treatment
- Validated against experiments

### RIXS Specialization:
- Momentum-resolved RIXS
- Polarization dependence
- dd and charge transfer excitations
- Magnetic excitations

### Modern Implementation:
- Python frontend
- Fortran backend
- MPI parallel
- Active development

## Inputs & Outputs
- **Input formats**:
  - Python scripts
  - Crystal field parameters
  - Slater integrals
  
- **Output data types**:
  - XAS spectra
  - RIXS maps
  - Eigenenergies and states
  - Transition matrix elements

## Installation
```bash
pip install edrixs
# Or with conda
conda install -c conda-forge edrixs
```

## Usage Examples
```python
import edrixs

# Set up atomic model for Ni2+ (d8)
case = edrixs.atom_solver.AtomSolver(
    shell_name='d',
    num_elec=8,
    slater_integrals=[...],
    spin_orbit_coupling=[...],
    crystal_field=[...]
)

# Calculate XAS spectrum
xas = case.get_xas(edge='L3')

# Calculate RIXS spectrum
rixs = case.get_rixs(
    incident_energy=852.0,
    pol_in=[1, 0, 0],
    pol_out=[0, 1, 0]
)
```

## Performance Characteristics
- **Speed**: Exact diagonalization scales with Hilbert space
- **Memory**: Limited by matrix size
- **Parallelization**: MPI for large calculations

## Limitations & Known Constraints
- **Model-based**: Requires crystal field parameters
- **Finite clusters**: Limited to small clusters
- **Scaling**: Exponential with number of orbitals
- **Parameter choice**: Requires expertise

## Comparison with Other Tools
- **vs CTM4XAS**: EDRIXS more flexible, open-source
- **vs Quanty**: Different implementation approaches
- **vs FEFF**: EDRIXS multiplet, FEFF real-space MS
- **Unique strength**: Open-source RIXS with Python API

## Application Areas
- Transition metal oxide RIXS
- Rare earth spectroscopy
- Magnetic excitations
- Crystal field analysis
- Charge transfer materials

## Best Practices
- Validate parameters against XAS first
- Use literature values for Slater integrals
- Check convergence with basis size
- Compare with experimental spectra

## Community and Support
- GitHub repository
- NSLS-II development
- CPC publication
- Active maintenance

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/NSLS-II/edrixs
2. Y.L. Wang et al., Comput. Phys. Commun. 243, 151 (2019)
3. Documentation: https://nsls-ii.github.io/edrixs/

**Confidence**: VERIFIED - Published in CPC

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GPL-3.0)
- Academic citations: Growing
- Active development: NSLS-II maintained
