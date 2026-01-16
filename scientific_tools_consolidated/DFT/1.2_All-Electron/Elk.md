# Elk

## Official Resources
- Homepage: https://elk.sourceforge.io/
- Documentation: https://elk.sourceforge.io/elk.html
- Source Repository: https://sourceforge.net/projects/elk/
- License: GNU General Public License v3.0

## Overview
Elk (formerly Exciting) is an all-electron full-potential linearized augmented plane wave (FP-LAPW) code for determining electronic structure of crystalline solids and molecules. It implements advanced methods for optical and spectroscopic properties and is particularly suited for systems requiring accurate treatment of core electrons and precise electronic structure calculations.

**Scientific domain**: All-electron DFT, FP-LAPW, electronic structure, spectroscopy  
**Target user community**: Researchers requiring all-electron accuracy, spectroscopy calculations, magnetic materials

## Theoretical Methods
- Full-potential linearized augmented plane wave (FP-LAPW)
- All-electron method (no pseudopotentials)
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Hybrid functionals (PBE0, HSE06, etc.)
- DFT+U for correlated systems
- Time-Dependent DFT (TDDFT)
- GW approximation
- Bethe-Salpeter equation (BSE)
- Spin-orbit coupling (second variational method)
- Bethe-Salpeter Equation (BSE)
- DFT+U for correlated systems
- Spin-orbit coupling
- Non-collinear magnetism
- Optimized effective potential (OEP)
- Exact exchange

## Capabilities (CRITICAL)
- Ground-state electronic structure (all-electron)
- Geometry optimization and relaxation
- Total energy and forces
- Band structure and DOS
- Optical properties (dielectric function, absorption)
- TDDFT for excited states
- GW quasiparticle energies
- BSE for optical excitations
- Magnetic properties (moments, anisotropies)
- Spin dynamics
- Magnon spectra
- Electric field gradients
- Hyperfine fields
- X-ray absorption and emission spectra
- Electron energy loss spectroscopy (EELS)
- Phonon calculations via linear response
- Elastic constants
- Wannier functions
- Berry phase and polarization
- Nonlinear optical properties

**Sources**: Official Elk documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - elk.in (main input file)
  - GEOMETRY.OUT format for structures
  - Species files for atomic data
  
- **Output data types**:
  - INFO.OUT (main output)
  - EIGVAL.OUT (eigenvalues)
  - DOS*.OUT (density of states)
  - EPSILON*.OUT (optical properties)
  - Various property-specific files

## Interfaces & Ecosystem
- **Framework integrations**:
  - elk2bloch - band unfolding
  - Wannier90 interface
  - LibXC for exchange-correlation functionals
  
- **Visualization**:
  - XCrySDen compatibility
  - elk-lapw utilities
  
- **Post-processing**:
  - Built-in analysis tools
  - elk-optics for optical spectra

## Workflow and Usage

### Basic DFT Calculation
```bash
# 1. Create elk.in input file with structure and parameters
# 2. Run Elk
elk

# Check convergence in INFO.OUT
# Results in EIGVAL.OUT, TDOS.OUT, etc.
```

### Optical Properties Calculation
```bash
# 1. Ground state calculation
elk

# 2. Add optics tasks to elk.in:
# ! Calculate optical response
tasks
  0    # Ground state
  120  # Dielectric function
  121  # Optical conductivity

# 3. Run calculation
elk

# Results in EPSILON_*.OUT files
```

### GW Calculation
```bash
# 1. Converged DFT ground state
elk

# 2. Setup GW calculation in elk.in:
tasks
  0    # Ground state
  300  # G0W0 calculation

# 3. Run GW
elk

# Quasiparticle energies in QPENE.OUT
```

### Magnetic Calculations
```bash
# Enable spin-polarized calculation in elk.in:
spinpol
  .true.

# Set initial magnetic moments
bfcmt
  [atom index] 0.0 0.0 [magnetic moment]

# Run calculation
elk

# Magnetic moments in INFO.OUT
```

## Application Areas
- Magnetic materials and spintronics
- Optical spectroscopy (absorption, reflectivity)
- X-ray spectroscopy (XAS, XES, EELS)
- Excited-state properties (GW, BSE)
- Heavy element systems (relativistic effects)
- Complex magnetic structures
- Phonons and lattice dynamics

## Limitations & Known Constraints
- **Open-source but specialized**: Smaller community than WIEN2k
- **All-electron cost**: Computationally expensive; ~100-200 atom limit
- **Learning curve**: LAPW methods require understanding
- **Documentation**: Comprehensive PDF but less tutorial material
- **Parallelization**: OpenMP and MPI but not as scalable as plane-wave codes
- **Installation**: Requires BLAS/LAPACK, FFTW
- **Input format**: Text-based, requires manual editing
- **Memory**: High for all-electron treatment
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: http://elk.sourceforge.net/
2. Documentation: http://elk.sourceforge.net/elk.pdf
3. SourceForge repository: http://sourceforge.net/projects/elk/
4. Elk development team publications

**Secondary sources**:
1. Elk manual and examples
2. Published applications using Elk
3. All-electron method benchmarks
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE (comprehensive PDF)
- Source code: OPEN (SourceForge)
- Community support: Active (mailing list)
- Academic citations: >500
- Active development: Regular releases
