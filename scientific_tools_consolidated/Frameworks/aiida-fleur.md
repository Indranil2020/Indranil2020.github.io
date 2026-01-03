# aiida-fleur

## Official Resources
- Homepage: https://github.com/JuDFTteam/aiida-fleur
- Documentation: https://aiida-fleur.readthedocs.io/
- Source Repository: https://github.com/JuDFTteam/aiida-fleur
- License: MIT License

## Overview
aiida-fleur is the AiiDA plugin for the FLEUR code, an all-electron DFT code based on the Full-Potential Linearized Augmented Plane Wave (FLAPW) method. It enables the automation of high-precision all-electron calculations, including magnetism and spin-orbit coupling, within the AiiDA framework.

**Scientific domain**: All-electron DFT, FLAPW, magnetism  
**Target user community**: FLEUR users, researchers requiring high-precision DFT

## Capabilities (CRITICAL)
- **Calculations**: Support for `fleur` and `inpgen` (input generator).
- **WorkChains**:
  - `FleurSCFWorkChain`: Robust self-consistency loop.
  - `FleurBandDosWorkChain`: Band structure and DOS.
  - `FleurRelaxWorkChain`: Geometry optimization.
  - `FleurMaeWorkChain`: Magnetic Anisotropy Energy.
- **Parsing**: Parses FLEUR XML output (`out.xml`) into AiiDA nodes.
- **Input Generation**: Pythonic interface to modify FLEUR parameters.

**Sources**: aiida-fleur documentation

## Inputs & Outputs
- **Input formats**: StructureData, Dict (parameters)
- **Output data types**: Dict (energies, charges), BandsData

## Interfaces & Ecosystem
- **FLEUR**: The backend engine
- **AiiDA**: Parent framework
- **JuDFT**: Developed by the Jülich DFT team

## Workflow and Usage
1. Define structure.
2. Submit `FleurSCFWorkChain`:
   ```python
   submit(WorkflowFactory('fleur.scf'), structure=structure, ...)
   ```
3. Automatic running of `inpgen` followed by `fleur` loops.

## Performance Characteristics
- Handles the complexity of all-electron convergence
- Specialized for magnetic materials

## Application Areas
- Magnetic materials (MAE, spin spirals)
- Oxides and f-electron systems (where all-electron is crucial)
- High-throughput screening of magnetic properties

## Community and Support
- Developed by Forschungszentrum Jülich (IAS-1)
- Active support

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/JuDFTteam/aiida-fleur
2. Documentation: https://aiida-fleur.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: All-electron DFT workflows, magnetism
