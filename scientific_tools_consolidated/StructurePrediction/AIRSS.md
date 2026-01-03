# AIRSS (Ab Initio Random Structure Searching)

## Official Resources
- Homepage: https://www.mtg.msm.cam.ac.uk/Codes/AIRSS
- Documentation: https://airss-docs.github.io/
- Source Repository: Distributed with license
- License: Academic License (Open source components available)

## Overview
AIRSS is a simple yet powerful method and software package for predicting crystal structures. The core idea is to generate random structures (sensibly constrained by density, symmetry, and distances) and relax them using ab-initio forces. Developed by Chris Pickard and collaborators, AIRSS is robust, easy to parallelize, and effective for a wide range of materials, especially under high pressure.

**Scientific domain**: Random structure searching, crystal prediction, materials discovery  
**Target user community**: Computational materials scientists, high-pressure physicists

## Theoretical Methods
- Random Structure Searching (RSS)
- Density functional theory relaxation
- Symmetry constraints (random space groups)
- Species swapping
- Shaking (for finding nearby minima)

## Capabilities (CRITICAL)
- Prediction of crystal structures from composition
- High-pressure phase prediction
- Point defect structure searching
- Surface and interface reconstruction
- Cluster structure search
- Constraint handling (distances, coordination)
- Integration with CASTEP (primary), VASP, Quantum ESPRESSO

**Sources**: AIRSS website, J. Phys.: Condens. Matter 23, 053201 (2011)

## Inputs & Outputs
- **Input formats**: `airss.pl` command line arguments, seed files
- **Output data types**: Relaxed structures (.res, .cell), energy rankings

## Interfaces & Ecosystem
- **CASTEP**: Deep integration and optimization
- **Command line**: Script-based workflow (`airss.pl`, `castep_relax`, etc.)
- **Analysis**: Tools for clustering and analyzing found structures

## Workflow and Usage
1. Generate random structures: `airss.pl -seed ...`
2. Relax structures: Submit to DFT code (e.g., CASTEP)
3. Collect results: Rank by enthalpy
4. Analyze: Identify unique structures and ground state

## Performance Characteristics
- Embarrassingly parallel (each random search is independent)
- High throughput
- Robust exploration of energy landscape
- Scales well to HPC clusters

## Application Areas
- High-pressure hydrides and oxides
- Defect structures in semiconductors
- Grain boundaries
- Battery materials
- Molecular crystals

## Community and Support
- Developed at University of Cambridge / UCL
- Active research use
- Workshops and tutorials

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.mtg.msm.cam.ac.uk/Codes/AIRSS
2. Documentation: https://airss-docs.github.io/
3. Publication: C. J. Pickard and R. J. Needs, J. Phys.: Condens. Matter 23, 053201 (2011)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Method: STANDARD (RSS)
- Development: ACTIVE (Pickard Group)
- Applications: Random structure searching, high pressure, defects, CASTEP integration
