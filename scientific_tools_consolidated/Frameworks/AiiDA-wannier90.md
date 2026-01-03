# AiiDA-wannier90

## Official Resources
- Homepage: https://github.com/aiidateam/aiida-wannier90
- Documentation: https://aiida-wannier90.readthedocs.io/
- Source Repository: https://github.com/aiidateam/aiida-wannier90
- License: MIT License

## Overview
AiiDA-wannier90 is the AiiDA plugin for Wannier90, the standard code for calculating Maximally Localized Wannier Functions (MLWFs). It enables automated wannierization workflows, handling the multi-step process (preprocessing, minimizing spread, post-processing) and facilitating the calculation of topological properties and transport.

**Scientific domain**: Wannier functions, electronic structure, topology  
**Target user community**: Users of Wannier90 and AiiDA

## Capabilities (CRITICAL)
- **Calculations**: Support for `wannier90.x` (main executable) and `postw90.x`.
- **WorkChains**: 
  - `Wannier90BaseWorkChain`: Runs Wannier90 with error handling.
  - `Wannier90BandsWorkChain`: Interpolated band structure calculation.
- **Automation**: Integration with SCDM (Selected Columns of the Density Matrix) method for automated initial projections.
- **Parsing**: Parses `.wout`, `.nnkp`, and other output files into AiiDA nodes.

**Sources**: AiiDA-wannier90 documentation

## Inputs & Outputs
- **Input formats**: AiiDA Nodes (StructureData, KpointsData, Dict parameters)
- **Output data types**: BandsData (interpolated bands), Dict (spreads, centers)

## Interfaces & Ecosystem
- **Wannier90**: The backend engine
- **AiiDA-QuantumESPRESSO**: Often used in conjunction (pw2wannier90 workflow)
- **AiiDA-VASP**: Can also be used with VASP inputs

## Workflow and Usage
1. Run DFT (e.g., QE).
2. Run `pw2wannier90` (if using QE).
3. Submit `Wannier90Calculation` or `Wannier90BandsWorkChain`.
   ```python
   submit(WorkflowFactory('wannier90.bands'), structure=structure, ...)
   ```

## Performance Characteristics
- Automates the complex wannierization pipeline
- Enables high-throughput wannierization via SCDM

## Application Areas
- Automated construction of tight-binding models
- High-throughput topological screening
- Berry phase calculations

## Community and Support
- Maintained by AiiDA and Wannier90 developers (EPFL, etc.)
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiidateam/aiida-wannier90
2. Documentation: https://aiida-wannier90.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Wannier90 workflows, automation
