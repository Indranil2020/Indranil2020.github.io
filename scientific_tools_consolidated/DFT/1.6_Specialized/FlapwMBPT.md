# FlapwMBPT (FLAPW with Many-Body Perturbation Theory)

## Official Resources
- **Principal Developer**: Andrey Kutepov (Brookhaven National Laboratory)
- **Host Group**: Comsuite (Rutgers/BNL correlation)
- **Source Repository**: Not publicly open-sourced; likely internal or available upon request (historical link `github.com/flapwmbpt/flapwmbpt` is 404).
- **Related Project**: [Comsuite](https://github.com/rutgersphysics/COMSUITE) (FlapwMBPT is often cited as the DFT/GW engine within this ecosystem).

## Overview
**FlapwMBPT** is an advanced electronic structure code developed by **Andrey Kutepov** (BNL) for over 30 years. It combines the rigorous **Full-Potential Linearized Augmented Plane Wave (FLAPW)** method with **Many-Body Perturbation Theory (MBPT)** to provide high-precision calculations of ground and excited states. It allows for consistent treatment of electronic correlations using GW approximation, self-consistent GW (scGW), and vertex corrections.

**Scientific domain**: Strongly Correlated Materials, Actinides, Spectroscopy, Methodological Development.
**Target user community**: Specialists in F-electron systems (Plutonium, Americium), GW method developers, and correlated electron researchers.

## Theoretical Methods
- **Basis Set**: FLAPW+LO (Full-Potential Linearized Augmented Plane Wave + Local Orbitals).
- **Relativity**: Scalar-relativistic and Fully Relativistic (Dirac) implementations.
- **Many-Body Theory**:
  - GW Approximation (G0W0, scGW).
  - Quasiparticle Self-Consistent GW (qsGW).
  - Vertex Corrections (GW+$\Gamma$).
  - Dynamic Mean Field Theory (DMFT) integration (via Comsuite).
- **Core Solvers**: All-electron treatment (no pseudopotentials).

## Capabilities
- **Electronic Structure**: Band structures, DOS, PDOS.
- **Excited States**: Quasiparticle energies, spectral functions.
- **Materials**: 
  - 3d Transition Metals.
  - Lanthanides and Actinides (f-electrons).
  - Transition Metal Oxides.
- **Thermodynamics**: Free energy calculations.
- **Parallelization**: MPI-based parallel execution.

## Key Strengths
- **All-Electron Accuracy**: Removes pseudopotential errors, critical for f-electrons and core states.
- **Self-Consistency**: Features one of the few fully self-consistent GW implementations.
- **Relativistic Effects**: High-accuracy treatment of spin-orbit coupling.
- **Vertex Corrections**: Advanced implementation beyond standard GW.

## Inputs & Outputs
- **Inputs**: Fortran-namelist style inputs (system definition, basis parameters, mixing schemes).
- **Outputs**:
  - Spectral functions ($A(\omega)$).
  - Band structures.
  - Self-energies ($\Sigma(\omega)$).
  - Thermodynamics data.

## Performance Characteristics
- **Scaling**: $N^3$ to $N^4$ depending on the level of theory (GW is expensive).
- **Hardware**: HPC-oriented, utilizing MPI for scaling across nodes.
- **Optimization**: Fortran90 codebase optimized for vector processors.

## Limitations & Known Constraints
- **Availability**: Unlike VASP or QE, it is not widely distributed; access is typically via collaboration or specific academic channels.
- **Cost**: All-electron GW is extremely computationally demanding compared to pseudopotential DFT.
- **Documentation**: Sparse public documentation; expertise required for usage.

## Comparison with Other Codes
- **vs WIEN2k**: Both are FLAPW. WIEN2k is the community standard for DFT; FlapwMBPT specializes in GW/MBPT.
- **vs exciting**: `exciting` is also all-electron GW (LAPW), but FlapwMBPT has distinct implementations of self-consistency and vertex corrections favored by the Kutepov/Gabriel Kotliar groups.
- **vs VASP/QE**: FlapwMBPT is all-electron (more accurate but slower) and focuses on many-body physics rather than high-throughput materials search.

## Verification & Sources
**Primary sources**:
1.  **BNL Profile**: [Andrey Kutepov Research](https://www.bnl.gov/) - Confirms development and ownership.
2.  **Comsuite**: Cited as a component code in the Comsuite/COMSCOPE project (Rutgers/BNL).
3.  **Literature**: Kutepov, A. L., et al. "Electronic structure of Pu and Am metals by self-consistent relativistic GW method." *Physical Review B* 85 (2012).

**Verification status**: ✅ VERIFIED (as Research Code)
- **Authenticity**: Confirmed existence and active research use.
- **Accessibility**: ⚠️ RESTRICTED (Public repo is missing/dead; code is likely shared privately).
