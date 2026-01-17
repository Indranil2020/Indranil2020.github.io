# Spex (Spectroscopy with EXcitations)

## Official Resources
- Homepage: https://www.flapw.de/spex/
- Documentation: https://www.flapw.de/spex/documentation/
- Source Repository: https://github.com/flapw-spex/spex
- License: Free for academic use

## Overview
Spex is an all-electron code for calculating quasiparticle energies and optical spectra using many-body perturbation theory (GW and Bethe-Salpeter equation). Developed primarily at Forschungszentrum Jülich, Spex uses the FLAPW method as input and implements sophisticated algorithms for GW self-energy and BSE kernel calculations. It is particularly powerful for accurate band gaps, quasiparticle band structures, and optical absorption spectra of solids.

**Scientific domain**: GW approximation, BSE, optical spectra, all-electron MBPT  
**Target user community**: Spectroscopy researchers, solid-state physicists, FLAPW users

## Theoretical Methods
- GW approximation (G₀W₀, GW₀, self-consistent GW)
- Bethe-Salpeter Equation (BSE)
- Random Phase Approximation (RPA)
- All-electron implementation
- FLAPW basis (uses FLEUR output)
- Full-frequency integration
- Contour deformation
- Plasmon-pole models

## Capabilities (CRITICAL)
- Quasiparticle energies (GW)
- Accurate band gaps
- Quasiparticle band structures
- Optical absorption spectra (BSE)
- Exciton binding energies
- Dielectric functions
- Photoemission spectra
- All-electron accuracy
- Core-level excitations
- Finite momentum transfer
- FLAPW input compatibility
- Production quality

**Sources**: Spex website (https://www.flapw.de/spex/)

## Key Strengths

### All-Electron GW:
- Full treatment of all electrons
- No pseudopotential approximations
- Core states included
- High accuracy
- Benchmark quality

### FLAPW Integration:
- Uses FLEUR output
- FLAPW basis advantages
- Full-potential accuracy
- All-electron wavefunctions
- Systematic approach

### Advanced Algorithms:
- Full-frequency GW
- Contour deformation
- Efficient RPA
- Optimized implementations
- Production performance

### Optical Spectra:
- BSE for excitons
- Accurate absorption
- Binding energies
- Oscillator strengths
- Experimental comparison

### Spectroscopy Focus:
- Photoemission (PES/IPES)
- Optical absorption
- Core-level excitations
- Finite momentum
- Comprehensive spectra

## Inputs & Outputs
- **Input formats**:
  - FLEUR wavefunctions and densities
  - Spex input files
  - k-point meshes
  - Frequency grids
  
- **Output data types**:
  - Quasiparticle energies
  - Band structures
  - Spectral functions
  - Optical spectra
  - Dielectric functions
  - BSE eigenstates

## Interfaces & Ecosystem
- **FLEUR Interface**:
  - Primary DFT input
  - FLAPW wavefunctions
  - Seamless integration
  - Tested workflow
  
- **Visualization**:
  - Standard plotting tools
  - Spectral data output
  - Band structure formats

## Workflow and Usage

### Typical Workflow:
1. Run FLEUR DFT calculation
2. Prepare Spex input
3. Run GW calculation
4. Analyze quasiparticle energies
5. Optional: Run BSE for optics
6. Extract and visualize spectra

### GW Calculation:
```bash
spex input.spex
# Computes GW corrections
```

### BSE for Optics:
- Calculate RPA dielectric function
- Solve BSE for excitons
- Obtain optical absorption spectrum

## Advanced Features

### GW Variants:
- G₀W₀ (one-shot)
- GW₀ (partially self-consistent)
- Self-consistent GW
- Different approximations
- User control

### Frequency Integration:
- Full-frequency approach
- Contour deformation
- Accurate self-energy
- No plasmon-pole approximation
- Systematic convergence

### BSE Implementation:
- Electron-hole interaction
- Exciton eigenstates
- Binding energies
- Oscillator strengths
- Finite momentum transfer

### All-Electron:
- Core electrons included
- Core-level excitations
- High-energy spectroscopy
- No frozen-core approximation
- Complete treatment

## Performance Characteristics
- **Speed**: Moderate (all-electron MBPT)
- **Accuracy**: Excellent (all-electron)
- **System size**: Unit cell to moderate
- **Scaling**: Standard GW scaling
- **Typical**: Research calculations

## Computational Cost
- **GW**: More expensive than DFT
- **BSE**: Additional cost for optics
- **All-electron**: Higher cost than pseudopotential
- **Accuracy**: Justifies computational expense
- **Production**: Feasible for research

## Limitations & Known Constraints
- **FLEUR dependency**: Requires FLEUR DFT input
- **System size**: Limited to moderate systems
- **Learning curve**: MBPT expertise needed
- **Computational cost**: All-electron expense
- **Platform**: Linux systems

## Comparison with Other Codes
- **vs BerkeleyGW**: Spex all-electron, BerkeleyGW pseudopotential
- **vs Yambo**: Both GW/BSE, Spex FLAPW-based
- **vs exciting**: Both all-electron, different algorithms
- **Unique strength**: FLAPW-based all-electron GW/BSE, full-frequency, FLEUR integration

## Application Areas

### Band Gap Corrections:
- Accurate fundamental gaps
- Quasiparticle bands
- Semiconductor properties
- Insulator gaps
- Band structure refinement

### Optical Spectroscopy:
- Absorption spectra
- Exciton physics
- Optical gaps
- Oscillator strengths
- Experimental comparison

### Photoemission:
- PES/IPES spectra
- Spectral functions
- Satellite features
- Comparison with experiments

### Materials Science:
- Electronic structure
- Excited states
- Optical properties
- Spectroscopy interpretation

## Best Practices

### DFT Preparation:
- Converged FLEUR calculation
- Appropriate k-mesh
- Sufficient empty states
- Quality wavefunctions

### GW Convergence:
- k-point convergence
- Frequency grid
- Empty states
- Cutoff parameters
- Systematic testing

### BSE Calculations:
- Appropriate transitions
- k-point sampling
- Exciton convergence
- Numerical parameters

## Community and Support
- Free for academic use
- FLAPW community
- Documentation available
- Research group support
- Publications and tutorials

## Educational Resources
- Spex documentation
- FLAPW school materials
- GW/BSE tutorials
- Published papers
- User examples

## Development
- Forschungszentrum Jülich
- Christoph Friedrich (main developer)
- Active development
- Regular updates
- Research-driven improvements

## Research Applications
- Accurate band gaps
- Optical spectra
- Quasiparticle physics
- Exciton studies
- Spectroscopy theory

## Technical Innovation

### All-Electron MBPT:
- No pseudopotentials
- Core states included
- Complete accuracy
- Benchmark calculations

### FLAPW Basis:
- Full-potential advantages
- Systematic basis
- All-electron treatment
- High precision

## FLAPW-GW Synergy
- FLEUR DFT input
- FLAPW wavefunctions
- All-electron consistency
- Integrated workflow
- Production quality

## Verification & Sources
**Primary sources**:
1. Spex website: https://www.flapw.de/spex/
2. GitHub: https://github.com/flapw-spex/spex
3. C. Friedrich et al., Comp. Phys. Comm. (2011)
4. Documentation and user manual

**Secondary sources**:
1. GW/BSE method literature
2. FLAPW method papers
3. Spectroscopy calculations
4. Application publications

**Confidence**: CONFIRMED - Established research code

**Verification status**: ✅ VERIFIED
- Website: ACCESSIBLE
- GitHub: Available
- Documentation: Comprehensive
- Community support: FLAPW/Jülich groups
- Active development: Regular updates
- Specialized strength: All-electron GW/BSE, FLAPW integration, full-frequency approach, optical spectra, photoemission, core excitations, benchmark accuracy
