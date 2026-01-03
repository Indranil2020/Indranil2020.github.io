# AFLOW (Automatic Flow)

## Official Resources
- Homepage: http://aflow.org/
- Documentation: http://aflow.org/documentation/
- Source Repository: http://aflow.org/install-aflow/ (Source available)
- License: GPL v3

## Overview
AFLOW (Automatic Flow) is a prominent high-throughput ab initio calculation framework and database. It automates the generation, simulation, and analysis of materials using DFT (primarily VASP). AFLOW manages a massive database of calculated properties (AFLOWlib) and provides tools for convex hull construction, prototype generation, and machine learning (AFLOW-ML).

**Scientific domain**: High-throughput DFT, materials database, materials informatics  
**Target user community**: Materials scientists, alloy designers

## Capabilities (CRITICAL)
- **Database**: AFLOWlib contains >3.5 million materials entries.
- **Automation**: Manages DFT calculations, error correction, and symmetry analysis.
- **Prototypes**: Extensive library of crystal prototypes for structure generation.
- **Properties**: Electronic structure, thermal properties (AEL/AGL for elastic/Debye), vibrations.
- **Machine Learning**: AFLOW-ML API for predicting properties (gap, bulk modulus, etc.) using trained models.
- **Symmtry**: AFLOW-SYM for symmetry analysis.

**Sources**: AFLOW website, Comp. Mater. Sci. 58, 218 (2012)

## Inputs & Outputs
- **Input formats**: `aflow.in` file, POSCAR
- **Output data types**: Database entries, properties JSON, web pages

## Interfaces & Ecosystem
- **VASP**: Primary DFT engine
- **Quantum ESPRESSO**: Support available
- **Web API**: REST API for querying the database programmatically

## Workflow and Usage
1. Generate `aflow.in` file from a structure or prototype.
2. Run `aflow --run` to execute the workflow (manages VASP).
3. Data is automatically post-processed and added to local or remote repository.
4. Query online database via API (`aflow --search`).

## Performance Characteristics
- Highly scalable (runs on top supercomputers)
- Massive database of pre-computed binary and ternary alloys

## Application Areas
- Alloy discovery (entropy stabilized alloys)
- Thermoelectric materials
- Superconductors
- Machine learning

## Community and Support
- Developed by Curtarolo Group (Duke University)
- Active development
- Regular schools/workshops

## Verification & Sources
**Primary sources**:
1. Homepage: http://aflow.org/
2. Publication: S. Curtarolo et al., Comp. Mater. Sci. 58, 218 (2012)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GPL)
- Development: ACTIVE
- Applications: High-throughput DFT, database, ML
