# AiiDA-QuantumESPRESSO

## Official Resources
- Homepage: https://github.com/aiidateam/aiida-quantumespresso
- Documentation: https://aiida-quantumespresso.readthedocs.io/
- Source Repository: https://github.com/aiidateam/aiida-quantumespresso
- License: MIT License

## Overview
AiiDA-QuantumESPRESSO is the official plugin for interfacing AiiDA with the Quantum ESPRESSO (QE) suite of codes (pw.x, ph.x, pp.x, etc.). It is one of the most mature and feature-rich plugins in the AiiDA ecosystem, providing robust workchains for standard DFT tasks like relaxations, band structures, and phonon calculations.

**Scientific domain**: DFT workflows, electronic structure, phonons  
**Target user community**: Quantum ESPRESSO users, AiiDA users

## Capabilities (CRITICAL)
- **Calculations**: Support for `pw.x` (SCF/Relax), `ph.x` (Phonons), `pp.x` (Post-processing), `dos.x`, `projwfc.x`, `cp.x` (Car-Parrinello), and more.
- **WorkChains**: 
  - `PwBaseWorkChain`: Robust handling of `pw.x` with automatic error recovery (restarting on convergence failure).
  - `PwRelaxWorkChain`: Geometry optimization workflow.
  - `PwBandsWorkChain`: Automated band structure calculation (seekpath integration).
  - `PhBaseWorkChain`: Phonon calculations.
- **Advanced Features**: Hubbard U handling, magnetism support.

**Sources**: AiiDA-QuantumESPRESSO documentation

## Inputs & Outputs
- **Input formats**: AiiDA Nodes (StructureData, KpointsData, UpfData, Dict)
- **Output data types**: BandsData, ArrayData (dos/projections), TrajectoryData, Dict (energies, forces)

## Interfaces & Ecosystem
- **AiiDA**: Core framework
- **Quantum ESPRESSO**: Calculation engine
- **aiida-pseudo**: Management of pseudopotentials (SSSP, PseudoDojo)

## Workflow and Usage
1. Load plugin: `CalculationFactory('quantumespresso.pw')`
2. Define inputs (structure, kpoints, parameters, pseudopotentials).
3. Submit WorkChain:
   ```python
   from aiida.plugins import WorkflowFactory
   PwRelax = WorkflowFactory('quantumespresso.pw.relax')
   submit(PwRelax, structure=structure, ...)
   ```

## Performance Characteristics
- Highly optimized for high-throughput
- "Protocol" system allows running standard calculations (moderate, precise, fast) with minimal input

## Application Areas
- High-throughput screening
- Phonon database generation
- Reproducible electronic structure data

## Community and Support
- Maintained by the AiiDA team (EPFL/PSI)
- Very active development and support

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiidateam/aiida-quantumespresso
2. Documentation: https://aiida-quantumespresso.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Quantum ESPRESSO workflows
