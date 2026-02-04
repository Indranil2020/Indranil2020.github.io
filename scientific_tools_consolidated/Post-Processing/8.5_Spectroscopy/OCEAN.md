# OCEAN (Obtaining Core Excitations from Ab initio electronic structure and NIST)

## Official Resources
- Homepage: https://www.nist.gov/services-resources/software/ocean
- Documentation: https://github.com/times-physics/ocean
- Source Repository: https://github.com/times-physics/ocean
- License: BSD 3-Clause License

## Overview
OCEAN is a code for calculating core-level spectra (XAS, XES, RIXS) from first principles using the Bethe-Salpeter Equation (BSE). It combines DFT (typically using Quantum ESPRESSO or ABINIT) with a GW/BSE approach to accurately treat core-hole interactions and excitonic effects. It is designed to handle periodic systems and provides accurate spectra for K-edges and L-edges.

**Scientific domain**: X-ray spectroscopy, BSE, core-level excitations  
**Target user community**: Spectroscopists, materials scientists

## Theoretical Methods
- Bethe-Salpeter Equation (BSE)
- Density Functional Theory (DFT)
- GW approximation (screening)
- Core-hole interaction
- Projector Augmented Wave (PAW) / Pseudopotentials
- Multiplet effects (limited)

## Capabilities (CRITICAL)
- Calculation of X-ray Absorption Spectra (XAS)
- X-ray Emission Spectra (XES)
- Resonant Inelastic X-ray Scattering (RIXS)
- Non-resonant Inelastic X-ray Scattering (NRIXS)
- Accurate treatment of core-hole screening
- Periodic systems (solids, surfaces)

**Sources**: OCEAN website, Comp. Phys. Comm. 182, 409 (2011)

## Key Strengths

### BSE Accuracy:
- Excitonic effects included
- Core-hole screening
- Many-body treatment
- Accurate near-edge

### First-Principles:
- No empirical parameters
- DFT-based workflow
- Periodic systems
- Multiple edges

### Open Source:
- BSD licensed
- GitHub hosted
- Active development
- NIST supported

## Inputs & Outputs
- **Input formats**: `ocean.in` (main input), DFT input files
- **Output data types**: Spectra files (energy vs intensity), absorption cross-sections

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Supported DFT backend
- **ABINIT**: Supported DFT backend
- **Parallelization**: MPI supported

## Workflow and Usage
1. Prepare DFT input structure.
2. Configure `ocean.in` (edges, screening parameters).
3. Run OCEAN script (automates DFT SCF, screening calculation, BSE Hamiltonian construction, and diagonalization).
4. Analyze spectral output.

## Performance Characteristics
- Computationally intensive (BSE solving)
- Scales with system size and basis set
- Parallelized for clusters

## Limitations & Known Constraints
- **Computational cost**: BSE is expensive
- **Convergence**: Requires careful k-point/band testing
- **Multiplet effects**: Limited treatment
- **Learning curve**: Complex setup

## Comparison with Other Tools
- **vs FEFF**: OCEAN BSE-based, FEFF real-space MS
- **vs exciting-XS**: Both BSE, different DFT backends
- **vs xspectra**: OCEAN more accurate, xspectra faster
- **Unique strength**: First-principles BSE for core levels

## Application Areas
- Transition metal oxides
- Battery materials (Li K-edge)
- Organic crystals
- Surface adsorbates

## Best Practices
- Converge k-points and bands
- Test screening parameters
- Validate against experiment
- Use appropriate core-hole treatment

## Community and Support
- Developed at NIST and University of Washington (Rehr group connection)
- Open-source (BSD)
- GitHub repository

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.nist.gov/services-resources/software/ocean
2. GitHub: https://github.com/times-physics/ocean
3. Publication: J. Vinson et al., Phys. Rev. B 83, 115106 (2011)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (NIST/Times Physics)
- Applications: BSE for X-ray spectroscopy
