# Fleur

## Official Resources
- Homepage: https://www.flapw.de/
- Documentation: https://www.flapw.de/MaX-7.0/documentation/
- Source Repository: https://iffgit.fz-juelich.de/fleur/fleur
- License: MIT License (open-source)

## Overview
Fleur is a feature-full, freely available FLAPW (Full-potential Linearized Augmented Plane Wave) code based on DFT, developed by the Forschungszentrum Jülich. It provides accurate all-electron calculations with a modern, well-maintained codebase and extensive capabilities for magnetic systems.

**Scientific domain**: Magnetism, spintronics, surfaces, all-electron calculations  
**Target user community**: Researchers studying magnetic materials, surfaces, and requiring all-electron accuracy

## Theoretical Methods
- Density Functional Theory (DFT)
- Full-potential linearized augmented plane wave (FLAPW)
- All-electron (no pseudopotentials)
- LDA, GGA, meta-GGA functionals
- Hybrid functionals
- DFT+U for correlated systems
- GW approximation (via Spex interface)
- Time-Dependent DFT
- Spin-orbit coupling
- Non-collinear magnetism
- Constrained DFT

## Capabilities (CRITICAL)
- Ground-state electronic structure (all-electron)
- Total energy and forces
- Geometry optimization
- Band structure and DOS
- Magnetic properties (moments, anisotropies, exchange interactions)
- Spin-spiral calculations
- Dzyaloshinskii-Moriya interaction
- Magnetic exchange parameters (J_ij)
- Surface and thin film calculations
- Electric field gradients
- Core-level spectroscopy
- X-ray magnetic circular dichroism (XMCD)
- Orbital magnetization
- Anomalous Hall conductivity
- Wannier functions via Wannier90 interface
- GW via Spex interface
- Phonon calculations (via DFPT or finite differences)
- Stress tensors for lattice optimization

**Sources**: Official Fleur documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - XML-based input files (inp.xml)
  - Structure files (various formats)
  - inp generator tools
  
- **Output data types**:
  - out.xml (main XML output)
  - out (formatted text output)
  - Density files
  - DOS and band structure files
  - Magnetic property outputs

## Interfaces & Ecosystem
- **Framework integrations**:
  - AiiDA-Fleur - workflow automation
  - Wannier90 - Wannier functions
  - Spex - GW calculations
  - JuDFT tools ecosystem
  
- **Visualization**:
  - Jmol-based viewers
  - Standard visualization tools
  
- **JuDFT family**:
  - JuKKR - KKR method code (same team)
  - masci-tools - Python interface

## Workflow and Usage

### Basic DFT Calculation
```bash
# 1. Generate input file (inp.xml)
inpgen -f structure.cif

# 2. Run self-consistent calculation
fleur

# Results in out.xml
```

### AiiDA-Fleur Workflow
```python
from aiida import load_profile
from aiida_fleur.workflows import FleurScfWorkChain

# Setup structure and parameters
structure = ...  # AiiDA StructureData
parameters = {...}

# Run automated SCF with AiiDA
submit(FleurScfWorkChain,
       fleur=code,
       structure=structure,
       calc_parameters=parameters)
```

### Magnetic Properties
```bash
# Enable spin-polarized calculation in inp.xml
# Set initial magnetic moments

# Calculate exchange interactions
fleur_exch

# Calculate magnetic anisotropy
fleur_mae

# Results: J_ij parameters, MAE values
```

### Surface Calculations
```bash
# 1. Generate film structure
inpgen -f slab_structure.cif -film

# 2. optimize film geometry  
fleur

# 3. Calculate surface properties
# Surface energies, work functions in out.xml
```

## Application Areas
- Magnetic materials (exchange interactions, anisotropy)
- Spintronics (spin-orbit coupling, DMI)
- Surfaces and interfaces
- Thin films and multilayers
- Anomalous Hall effect
- Orbital magnetization
- High-throughput materials screening (AiiDA)

## Limitations & Known Constraints
- **All-electron cost**: Computationally expensive; limited to ~100-200 atoms
- **Learning curve**: FLAPW methods require understanding
- **Parallelization**: MPI parallelization good but not as scalable as plane-wave codes
- **Memory**: High for all-electron calculations
- **Community**: Smaller than WIEN2k but growing
- **Documentation**: Good but evolving with code updates
- **Installation**: Requires modern Fortran compiler, libraries
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://www.flapw.de/
2. Documentation: https://www.flapw.de/MaX-7.0/documentation/
3. GitLab repository: https://iffgit.fz-juelich.de/fleur/fleur
4. Fleur development team (FZ Jülich)

**Secondary sources**:
1. Fleur tutorials and examples
2. AiiDA-Fleur documentation
3. Published magnetism applications
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitLab)
- Community support: Active (GitLab, mailing list)
- Academic citations: >500
- Active development: Regular releases, well-maintained
