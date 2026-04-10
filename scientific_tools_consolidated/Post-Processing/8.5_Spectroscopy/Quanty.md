# Quanty

## Official Resources
- Homepage: https://www.quanty.org/
- Documentation: https://www.quanty.org/documentation/
- Source Repository: https://www.quanty.org/download/
- License: Free for academic use

## Overview
**Quanty** is a many-body script language based on Lua designed for the calculation of X-ray spectroscopy including XAS, XES, RIXS, NIXS, and XPS. It uses exact diagonalization and Lanczos methods to compute multiplet spectra for correlated materials, ranging from crystal field theory to ligand field theory and post-DFT treatments.

**Scientific domain**: X-ray spectroscopy, correlated electron systems, multiplet calculations  
**Target user community**: Researchers at synchrotron facilities studying transition metal and rare-earth X-ray spectra

## Theoretical Methods
- Exact diagonalization (ED) of many-body Hamiltonians
- Lanczos algorithm for large Hilbert spaces
- Crystal field theory
- Ligand field theory (charge transfer models)
- Configuration interaction (CI)
- Post-DFT treatment of core-level spectroscopy
- Spin-orbit coupling
- Multiplet theory (Cowan's code integration)
- Determinant-based spectroscopy

## Capabilities (CRITICAL)
- X-ray absorption spectroscopy (XAS)
- X-ray emission spectroscopy (XES)
- Resonant inelastic X-ray scattering (RIXS)
- Non-resonant inelastic X-ray scattering (NIXS)
- X-ray photoelectron spectroscopy (XPS)
- X-ray magnetic circular dichroism (XMCD)
- RIXS-MCD (magnetic circular dichroism in RIXS)
- Partial fluorescence yield (PFY)
- Inverse partial fluorescence yield (IPFY)
- Ground-state properties (density matrices)
- Temperature-dependent spectra
- Polarization-dependent spectra

**Sources**: Official Quanty documentation, J. Phys. Conf. Ser. 712, 012005 (2016)

## Key Strengths

### Many-Body Treatment:
- Full multiplet calculations
- Beyond independent-particle approximation
- Proper treatment of core-hole effects
- Spin-orbit coupling included
- Handles strong correlation

### Flexible Scripting:
- Lua-based input language
- User-defined Hamiltonians
- Custom basis sets
- Modular calculation flow
- Extensible for new spectroscopies

### Broad Spectroscopy Coverage:
- XAS, XES, RIXS, NIXS, XPS
- XMCD and RIXS-MCD
- Multiple edge types (K, L, M)
- Polarization dependence
- Temperature dependence

### DFT Integration:
- Interfaces with DFT for radial integrals
- Crystal field parameters from DFT
- Slater-Condon parameters from DFT
- Charge transfer energies from DFT
- Combined DFT+multiplet approach

## Inputs & Outputs
- **Input formats**:
  - Lua script files
  - Hamiltonian definitions
  - Basis set specifications
  - Spectroscopy parameters
  
- **Output data types**:
  - Spectra (XAS, XES, RIXS cross-sections)
  - Density matrices
  - Ground state energies
  - Transition matrix elements
  - Dichroism spectra

## Interfaces & Ecosystem
- **Crispy**: GUI for Quanty (https://github.com/mretegan/crispy)
- **Quanty4RIXS**: MATLAB interface for RIXS calculations
- **DFT codes**: Parameter extraction from DFT
- **Cowan's code**: Atomic multiplet parameters

## Performance Characteristics
- **Speed**: Lanczos algorithm enables large Hilbert spaces
- **Accuracy**: Full multiplet treatment, beyond single-particle
- **System size**: Limited by Hilbert space dimension (typically d/f electron systems)
- **Parallelization**: Limited native parallelization

## Computational Cost
- **XAS**: Moderate (Lanczos)
- **RIXS**: Expensive (two-particle Green's function)
- **Full multiplet**: Scales with Hilbert space dimension
- **Typical**: Minutes to hours per spectrum

## Limitations & Known Constraints
- **Hilbert space**: Limited by exact diagonalization
- **Basis size**: Restricts to few correlated sites
- **No periodic boundary conditions**: Cluster/impurity models only
- **Learning curve**: Lua scripting required
- **Radial integrals**: Need external DFT or Cowan's code
- **No self-consistency**: Parameters from external sources

## Comparison with Other Codes
- **vs CTM4XAS**: Quanty is more flexible, CTM4XAS is simpler GUI
- **vs FEFF**: Quanty handles multiplet effects, FEFF is single-particle
- **vs EDRIXS**: Both use ED; Quanty is Lua-based, EDRIXS is Python
- **vs Crispy**: Crispy is a GUI wrapper around Quanty
- **Unique strength**: Flexible many-body scripting for X-ray spectroscopy, broad spectroscopy coverage, Lua-based extensibility

## Application Areas

### Transition Metal Spectroscopy:
- 3d transition metal L-edge XAS/RIXS
- Crystal field and charge transfer effects
- Spin-state transitions
- Mixed-valence systems

### Rare-Earth Spectroscopy:
- 4f rare-earth M-edge spectroscopy
- Multiplet structures
- Crystal field splitting
- Magnetic properties

### Strongly Correlated Materials:
- Cuprates and manganites
- Charge transfer insulators
- Mott insulators
- High-Tc superconductors

### Catalysts and Functional Materials:
- Active site characterization
- Oxidation state determination
- Spin-state analysis
- Coordination environment

## Best Practices

### Parameter Selection:
- Use DFT for crystal field parameters
- Calibrate Slater-Condon reductions
- Validate against experimental spectra
- Consider charge transfer models

### RIXS Calculations:
- Start with XAS to validate parameters
- Use appropriate k-weights
- Monitor convergence of Lanczos
- Consider temperature effects

### Hamiltonian Construction:
- Include relevant orbitals only
- Use symmetry to reduce Hilbert space
- Validate crystal field parameters
- Check spin-orbit coupling strength

## Community and Support
- Free for academic use
- Active development
- Documentation and tutorials on quanty.org
- Used at major synchrotron facilities (ESRF, ALS, NSLS-II)
- Quanty4RIXS MATLAB companion

## Verification & Sources
**Primary sources**:
1. Official website: https://www.quanty.org/
2. M. W. Haverkort, J. Phys. Conf. Ser. 712, 012005 (2016)
3. Crispy GUI: https://github.com/mretegan/crispy
4. Quanty4RIXS: J. Synchrotron Rad. 25, 904 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: Available for download
- Community support: Active (synchrotron community)
- Academic citations: >100 (main papers)
- Active development: Ongoing
- Specialized strength: Many-body X-ray spectroscopy, Lua scripting, full multiplet treatment
