# PyLada

## Official Resources
- Homepage: https://github.com/pylada/pylada-light
- Documentation: Minimal
- Source Repository: https://github.com/pylada/pylada-light
- License: BSD License

## Overview
PyLada is a Python framework for high-throughput computational materials science. It provides tools for managing DFT calculations (wrapping VASP and Crystal), handling crystal structures, and managing job submission to clusters. It is designed to be lightweight and scriptable.

**Scientific domain**: Workflow management, high-throughput DFT  
**Target user community**: VASP/Crystal users

## Capabilities (CRITICAL)
- **Wrappers**: Interfaces for VASP and Crystal.
- **Structure**: Crystal structure classes (superlattices, etc.).
- **Job Management**: Simple job folder management.
- **High-throughput**: Tools for enumerating structures (e.g., alloys).

**Sources**: PyLada GitHub

## Inputs & Outputs
- **Input formats**: Python scripts
- **Output data types**: VASP outputs

## Interfaces & Ecosystem
- **VASP**: Primary engine.
- **Crystal**: Supported.

## Workflow and Usage
1. Define structure.
2. Setup functional: `vasp = Vasp()`
3. Launch job.

## Performance Characteristics
- Lightweight Python wrapper.

## Application Areas
- Alloy theory (Cluster expansion)
- Defect calculations

## Community and Support
- Developed at NREL (Lany group)
- **Status**: Less active than Pymatgen/Atomate.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/pylada/pylada-light

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: MINIMAL
- Source: OPEN (GitHub)
- Development: SLOW/MAINTENANCE
- Applications: VASP wrapper
