# QMCPACK-addons (Nexus)

## Official Resources
- Homepage: https://qmcpack.org/
- Documentation: https://nexus-workflow.readthedocs.io/
- Source Repository: https://github.com/QMCPACK/qmcpack (nexus directory)
- License: UIUC/NCSA Open Source License

## Overview
Nexus is a workflow management system distributed with QMCPACK. It automates the complex series of steps required for Quantum Monte Carlo (QMC) calculations: generating trial wavefunctions (from DFT codes like QE, VASP), optimizing Jastrow factors, and running DMC. It manages dependencies and job submission.

**Scientific domain**: QMC workflows  
**Target user community**: QMCPACK users

## Capabilities (CRITICAL)
- **Automation**: DFT -> Conversion -> Opt -> DMC pipeline.
- **Codes**: Supports QMCPACK, Quantum ESPRESSO, VASP, PySCF, GAMESS.
- **System**: Structure generation (supercells).
- **Analysis**: Statistical analysis of QMC data.

**Sources**: Nexus documentation, Comp. Phys. Comm. 203, 110 (2016)

## Inputs & Outputs
- **Input formats**: Python script
- **Output data types**: QMCPACK inputs, job scripts

## Interfaces & Ecosystem
- **QMCPACK**: Primary target.
- **DFT Codes**: QE, VASP, etc.

## Workflow and Usage
1. Write Python script defining physical system and simulation cascade.
2. `run_project()`

## Performance Characteristics
- Simplifies the tedious setup of QMC.

## Application Areas
- High-accuracy electronic structure.
- Reference data generation.

## Community and Support
- Developed by Krogel (ORNL) and QMCPACK team.

## Verification & Sources
**Primary sources**:
1. Documentation: https://nexus-workflow.readthedocs.io/
2. Publication: J. T. Krogel, Comp. Phys. Comm. 203, 110 (2016)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN
- Development: ACTIVE
- Applications: QMC workflows
