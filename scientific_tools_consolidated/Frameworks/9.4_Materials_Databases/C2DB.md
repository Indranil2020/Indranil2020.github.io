# C2DB (Computational 2D Materials Database)

## Official Resources
- Homepage: https://cmr.fysik.dtu.dk/c2db/c2db.html
- Documentation: https://cmr.fysik.dtu.dk/c2db/c2db.html
- Source Repository: Data available via ASE/ASR
- License: Open Data

## Overview
C2DB is a highly curated database of two-dimensional materials calculated using the Atomic Simulation Environment (ASE) and the GPAW electronic structure code. It systematically classifies and characterizes thousands of 2D monolayers, calculating a wide range of properties including stability, stiffness, topological invariants, and excited state properties (GW band gaps, excitons).

**Scientific domain**: 2D materials, high-throughput DFT, many-body perturbation theory  
**Target user community**: 2D materials researchers

## Capabilities (CRITICAL)
- **Database**: >4000 2D materials.
- **Stability**: Dynamical stability (phonons) and thermodynamic stability (convex hull).
- **Electronic**: PBE and HSE06 band structures, GW quasiparticle gaps.
- **Optical**: BSE absorption spectra (excitons).
- **Topology**: Berry phases, Chern numbers, Z2 invariants.
- **Magnetic**: Heisenberg exchange parameters, magnetic anisotropy.
- **Mechanical**: Elastic stiffness tensors.

**Sources**: C2DB website, 2D Mater. 8, 025021 (2021)

## Inputs & Outputs
- **Input formats**: Searchable web interface
- **Output data types**: ASE database format (.db), JSON

## Interfaces & Ecosystem
- **ASR (Atomic Simulation Recipes)**: The workflow engine used to generate C2DB.
- **GPAW**: The primary calculation code.
- **ASE**: The underlying framework.
- **MyQueue**: Used for job scheduling.

## Workflow and Usage
1. **Web**: Browse materials at cmr.fysik.dtu.dk.
2. **Python**: Download the database via `ase`:
   ```python
   from ase.db import connect
   db = connect('c2db.db')
   ```

## Performance Characteristics
- High-fidelity data (GW/BSE)
- Focus exclusively on 2D materials

## Application Areas
- Discovery of new 2D semiconductors
- Topological insulators
- Single-photon emitters
- Catalysis

## Community and Support
- Developed by CAMD (Thygesen Group) at DTU
- Open access

## Verification & Sources
**Primary sources**:
1. Homepage: https://cmr.fysik.dtu.dk/c2db/c2db.html
2. Publication: S. Haastrup et al., 2D Mater. 5, 042002 (2018); M. N. Gjerding et al., 2D Mater. 8, 025021 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (Recipes via ASR)
- Development: ACTIVE
- Applications: 2D materials database, GW/BSE
