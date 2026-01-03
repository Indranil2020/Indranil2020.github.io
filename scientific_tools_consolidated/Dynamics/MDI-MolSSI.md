# MDI drivers (MolSSI Driver Interface)

## Official Resources
- Homepage: https://molssi-mdi.github.io/MDI_Library/
- Documentation: https://molssi-mdi.github.io/MDI_Library/
- Source Repository: https://github.com/MolSSI-MDI/MDI_Library
- License: BSD 3-Clause License

## Overview
The MolSSI Driver Interface (MDI) is a standardized API that allows different computational chemistry codes to communicate and exchange data during runtime. It enables interoperability between codes (e.g., a quantum chemistry code and a molecular dynamics driver) without requiring them to be linked into a single executable. MDI drivers facilitate complex workflows like QM/MM, advanced sampling, and machine learning integration.

**Scientific domain**: Interoperability, code coupling, multiscale modeling, QM/MM  
**Target user community**: Developers of computational chemistry software, researchers needing coupled codes

## Theoretical Methods
- Client-Server Architecture
- Runtime Data Exchange
- Coupled Simulations (QM/MM, AIMD)
- Force Bridging
- External Driver Control

## Capabilities (CRITICAL)
- Standardized communication between simulation codes
- Language-agnostic interface (C, C++, Fortran, Python)
- Support for MPI parallelization within and between codes
- Enables QM/MM without monolithic codebases
- Allows Python drivers to control compiled HPC codes
- Dynamic library loading or TCP/IP socket communication

**Sources**: MDI documentation, J. Chem. Theory Comput. (MolSSI publications)

## Inputs & Outputs
- **Input formats**: MDI command-line options (`-mdi`), driver scripts
- **Output data types**: Exchange of coordinates, forces, energies, virials, charges, etc.

## Interfaces & Ecosystem
- **Supported Codes (MDI-compliant)**: LAMMPS, Quantum ESPRESSO, Psi4, Q-Chem, Molcas, OpenMM, Tinker, etc.
- **Python**: MDI Python package for writing drivers
- **Standards**: Defined by MolSSI (Molecular Sciences Software Institute)

## Workflow and Usage
1. Launch code A (e.g., QM code) in MDI driver mode: `psi4 -mdi "role=DRIVER ..."`
2. Launch code B (e.g., MD code) in MDI engine mode: `lmp_mpi -mdi "role=ENGINE ..."`
3. Codes handshake and exchange data as defined by the MDI standard
4. Driver controls the simulation loop

## Performance Characteristics
- Minimal overhead for communication
- Enables heterogeneous computing (e.g., GPU MD + CPU QM)
- Flexible coupling strategies (tight vs. loose)

## Application Areas
- QM/MM simulations
- Ab initio molecular dynamics (AIMD)
- Machine learning potentials coupled to MD
- Advanced sampling (Metadynamics driven by external code)
- Multiscale modeling

## Community and Support
- Developed by MolSSI (US NSF-funded institute)
- Active development and standardization
- Workshops and hackathons
- Growing adoption in major codes

## Verification & Sources
**Primary sources**:
1. Homepage: https://molssi-mdi.github.io/MDI_Library/
2. GitHub: https://github.com/MolSSI-MDI/MDI_Library

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (MolSSI)
- Applications: Interoperability, code coupling, standard API
