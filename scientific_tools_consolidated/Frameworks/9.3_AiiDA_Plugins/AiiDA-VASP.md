# AiiDA-VASP

## Official Resources
- Homepage: https://github.com/aiida-vasp/aiida-vasp
- Documentation: https://aiida-vasp.readthedocs.io/
- Source Repository: https://github.com/aiida-vasp/aiida-vasp
- License: MIT License

## Overview
AiiDA-VASP is the official AiiDA plugin for the Vienna Ab initio Simulation Package (VASP). It provides an interface to run VASP calculations within the AiiDA workflow engine. It handles input generation (INCAR, POSCAR, KPOINTS, POTCAR), job submission, and parsing of output files (XML, OUTCAR) into the AiiDA database, preserving full provenance.

**Scientific domain**: DFT workflows, high-throughput VASP  
**Target user community**: VASP users using AiiDA

## Capabilities (CRITICAL)
- **Calculations**: Support for standard VASP runs (`VaspCalculation`).
- **WorkChains**: High-level workflows for relaxation (`RelaxWorkChain`), bands (`BandsWorkChain`), and convergence.
- **Parsing**: Robust parsing of VASP outputs into Python dictionaries and AiiDA nodes.
- **POTCAR Management**: Tools for managing VASP pseudopotentials families.
- **Error Handling**: Basic error handling via AiiDA mechanisms.

**Sources**: AiiDA-VASP documentation

## Inputs & Outputs
- **Input formats**: AiiDA Nodes (StructureData, KpointsData, Dict for parameters)
- **Output data types**: BandsData, TrajectoryData, Dict (energies, forces)

## Interfaces & Ecosystem
- **AiiDA**: The parent framework
- **VASP**: The calculation engine
- **Pymatgen**: Often used for underlying IO logic

## Workflow and Usage
1. Install plugin: `pip install aiida-vasp`
2. Configure code and potentials.
3. Submit WorkChain:
   ```python
   builder = WorkflowFactory('vasp.relax').get_builder()
   builder.structure = structure
   builder.parameters = Dict(dict={'incar': {...}})
   submit(builder)
   ```

## Performance Characteristics
- Dependent on VASP performance
- AiiDA overhead is minimal compared to DFT runtime

## Application Areas
- High-throughput screening of materials properties
- Database generation
- reproducible VASP calculations

## Community and Support
- Maintained by AiiDA developers and community contributors
- Active GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiida-vasp/aiida-vasp
2. Documentation: https://aiida-vasp.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: VASP interface for AiiDA
