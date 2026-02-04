# exciting-XS

## Official Resources
- Homepage: https://exciting-code.org/
- Documentation: https://exciting-code.org/carbon-xs-spectroscopy
- Source Repository: https://github.com/exciting-code/exciting
- License: GNU General Public License v2.0

## Overview
exciting-XS refers to the X-ray spectroscopy capabilities within the **exciting** code. It is a full-potential all-electron DFT code based on the linearized augmented plane-wave (LAPW) method. It can calculate core-level spectra (XAS, XES, EELS) using the Bethe-Salpeter Equation (BSE) or Time-Dependent DFT (TDDFT), providing highly accurate descriptions of core excitations including excitonic effects.

**Scientific domain**: X-ray spectroscopy, core-level excitations, all-electron DFT  
**Target user community**: Spectroscopists, materials scientists requiring high precision

## Theoretical Methods
- Full-potential LAPW+lo method
- Bethe-Salpeter Equation (BSE) for core states
- Time-Dependent Density Functional Theory (TDDFT)
- Core-hole interaction (electron-hole attraction)
- Spin-orbit coupling in core states
- All-electron treatment (no pseudopotentials)

## Capabilities (CRITICAL)
- Calculation of K, L, M edge XAS and XES
- Solution of the BSE for core excitations (accurate excitonic effects)
- Electron Energy Loss Spectroscopy (EELS) at core edges
- Momentum-dependent excitations (q-dependence)
- Real-space visualization of core excitons
- High-precision results due to LAPW basis

**Sources**: exciting documentation, Comp. Phys. Comm. 185, 2080 (2014)

## Key Strengths

### All-Electron Accuracy:
- Full-potential LAPW
- No pseudopotential errors
- Core state treatment
- High precision

### BSE for Core Levels:
- Excitonic effects
- Core-hole interaction
- Momentum dependence
- Exciton visualization

### Open Source:
- GPL licensed
- GitHub hosted
- Active development
- Good documentation

## Inputs & Outputs
- **Input formats**: `input.xml` (comprehensive XML configuration)
- **Output data types**: `XAS_*.xml` (spectra), `EPSILON_*.OUT` (dielectric function), `EXCITON_*.OUT`

## Interfaces & Ecosystem
- **Visualization**: Tools to plot XML output
- **Cluster**: Parallelized (MPI/OpenMP) for HPC
- **Tools**: `excitingscripts` for post-processing

## Workflow and Usage
1. Perform ground state SCF calculation.
2. Configure `input.xml` for `xs` calculation (define edge, screening, BSE parameters).
3. Run `exciting_smp` or `exciting_serial`.
4. Analyze spectral files.

## Performance Characteristics
- Computationally demanding (BSE scaling)
- Requires convergence of empty states
- Highly accurate but slower than pseudopotential methods

## Limitations & Known Constraints
- **Computational cost**: BSE is expensive
- **Empty states**: Requires convergence
- **LAPW complexity**: Steeper learning curve
- **System size**: Limited by cost

## Comparison with Other Tools
- **vs OCEAN**: Both BSE, exciting all-electron
- **vs FEFF**: exciting periodic, FEFF cluster
- **vs xspectra**: exciting more accurate, xspectra faster
- **Unique strength**: All-electron BSE for core levels

## Application Areas
- Core-level spectroscopy of complex oxides
- Excitonic effects in insulators and semiconductors
- Anisotropic X-ray absorption
- Core-hole screening physics

## Best Practices
- Converge empty states carefully
- Use appropriate k-point mesh
- Test BSE parameters
- Validate with experiment

## Community and Support
- Open-source (GPL)
- Developed by exciting team (Humboldt Univ. Berlin, etc.)
- Active forum and tutorials

## Verification & Sources
**Primary sources**:
1. Homepage: https://exciting-code.org/
2. Publication: A. Gulans et al., J. Phys.: Condens. Matter 26, 363202 (2014)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: BSE for core levels, all-electron XAS
