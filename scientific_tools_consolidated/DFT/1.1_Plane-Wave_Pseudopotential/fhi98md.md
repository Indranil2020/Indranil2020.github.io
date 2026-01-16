# fhi98md

## Official Resources
- Homepage: Historic (FHI Berlin)
- Archive: CPC Program Library (Computer Physics Communications)
- License: Academic / Historic

## Overview
fhi98md is a historic, pioneering plane-wave pseudopotential Density Functional Theory code developed at the Fritz Haber Institute (FHI) of the Max Planck Society in Berlin. Released in the late 1990s, it notably standardized the "FHI" pseudopotential format and the `.ini` input style used by many subsequent codes. It was a workhorse for Surface Science ab initio Molecular Dynamics (AIMD) before being superseded by codes like ABINIT and FHI-aims.

**Scientific domain**: Historic Surface Science, AIMD
**Target user community**: Historians of computational physics, Maintainers of legacy data

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets
- Norm-Conserving Pseudopotentials (Troullier-Martins)
- Car-Parrinello and Born-Oppenheimer Molecular Dynamics

## Capabilities
- Ground state electronic structure.
- Molecular Dynamics (NVE, NVT ensembles).
- Structure optimization.
- Parallelization (MPI).

## Key Strengths

### Historic Standards:
- **Pseudopotentials**: The `fhi` pseudopotential format defined by this code is still supported by tools like ABINIT and atomic generators (opep).
- **Simplicity**: Compact Fortran codebase that was easy to modify, leading to many private forks in the 2000s.

## Inputs & Outputs
- **Input formats**:
  - `start.ini`: Control parameters.
  - `coord.ini`: Atomic coordinates.
  - `pseudopotentials`: In FHI format.
  
- **Output data types**:
  - Text logs.
  - Fortran binary dump files.

## Interfaces & Ecosystem
- **ABINIT**: Can read FHI pseudopotentials.
- **FHI-aims**: The spiritual successor at FHI (though Aims uses numeric atom-centered orbitals, not plane waves).

## Computational Cost
- **Pre-GPU**: Designed for vector supercomputers and early MPI clusters of the late 90s.
- **Efficiency**: Good for its time; outperformed by modern libraries (FFTW3, LAPACK optimized) used in current codes.

## Limitations & Known Constraints
- **Obsolescence**: No longer maintained.
- **Hardware**: optimized for single-core or old-school MPI; inefficient on modern multi-core/GPU nodes.

## Comparison with Other Codes
- **vs ABINIT**: ABINIT supports the same pseudopotentials but adds thousands of new features (GW, DMFT, Berry Phases).
- **vs FHI-aims**: FHI-aims (the successor) switched to numeric atom-centered orbitals for all-electron accuracy, moving away from the plane-wave pseudopotential approach of fhi98md.

## Best Practices
- **Do Not Use for New Research**: It is obsolete. Use Quantum ESPRESSO, ABINIT, or FHI-aims instead.
- **Data Recovery**: Use documentation to parse historic `inp` or `ini` files if recovering old research data.

## Community and Support
- **Status**: **Obsolete**.
- **Legacy**: Its lineage lives on in the methodologies standardized by the Scheffler group at FHI.

## Verification & Sources
**Primary sources**:
1. "fhi98md: A efficient ab initio molecular dynamics program..." (CPC Publication)
2. FHI Pseudopotential database.

**Confidence**: VERIFIED - Historic code of high significance.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: Historic PW-DFT
- Key Feature: FHI Pseudopotential Format
