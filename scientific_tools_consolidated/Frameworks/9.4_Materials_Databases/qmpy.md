# qmpy (Quantum Materials Python)

## Official Resources
- Homepage: https://github.com/wolverton-research-group/qmpy
- Documentation: https://qmpy.readthedocs.io/
- Source Repository: https://github.com/wolverton-research-group/qmpy
- License: MIT License

## Overview
qmpy is the backend software stack for the Open Quantum Materials Database (OQMD). It provides tools for creating, managing, and analyzing high-throughput DFT calculations. It includes modules for structure manipulation, VASP input generation, job management, and thermodynamic analysis (phase diagrams).

**Scientific domain**: High-throughput DFT, thermodynamics, OQMD backend  
**Target user community**: OQMD users, Wolverton group, alloy researchers

## Capabilities (CRITICAL)
- **OQMD Interface**: Direct interaction with OQMD database models (Django-based).
- **Phase Diagrams**: Robust construction of convex hulls and stability analysis.
- **Workflow**: Manages VASP jobs (local or cluster).
- **Analysis**: Formation energies, reaction energies, grand canonical stability.

**Sources**: qmpy documentation

## Inputs & Outputs
- **Input formats**: Structure files, database queries
- **Output data types**: Phase diagrams, VASP inputs

## Interfaces & Ecosystem
- **Django**: Uses Django ORM for database abstraction.
- **VASP**: Primary calculation engine.
- **OQMD**: The data it manages.

## Workflow and Usage
1. Setup database settings.
2. Ingest structure: `entry = Entry.create("POSCAR")`
3. Analyze stability: `phase_diagram = PhaseDiagram(entry.elements); phase_diagram.stability(entry)`

## Performance Characteristics
- Database-centric (MySQL/PostgreSQL)
- Optimized for thermodynamic analysis

## Application Areas
- Maintaining OQMD
- Local high-throughput studies
- Phase stability research

## Community and Support
- Developed by Wolverton Group (Northwestern)
- Open source

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/wolverton-research-group/qmpy

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: OQMD, thermodynamics, VASP workflow
