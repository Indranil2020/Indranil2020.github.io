# AiiDA-yambo

## Official Resources
- Homepage: https://github.com/yambo-code/aiida-yambo
- Documentation: https://aiida-yambo.readthedocs.io/
- Source Repository: https://github.com/yambo-code/aiida-yambo
- License: MIT License

## Overview
AiiDA-yambo is the AiiDA plugin for the Yambo code, which performs Many-Body Perturbation Theory (MBPT) calculations (GW approximation, Bethe-Salpeter Equation). It allows users to automate complex excited-state calculations, including convergence tests for the many numerical parameters involved in GW/BSE.

**Scientific domain**: Many-body perturbation theory, excited states, GW, BSE  
**Target user community**: Yambo users, AiiDA users

## Capabilities (CRITICAL)
- **Calculations**: Support for `yambo` (initialization), `yambo` (run), `p2y` (interface).
- **WorkChains**:
  - `YamboConvergence`: Automated convergence of parameters (k-points, bands, cutoff, FFT grid).
  - `YamboRestart`: Automatic handling of walltime and errors.
  - `YamboWorkflow`: End-to-end GW calculation starting from DFT.
- **Parsing**: Extraction of quasiparticle energies and optical spectra.

**Sources**: AiiDA-yambo documentation

## Inputs & Outputs
- **Input formats**: AiiDA Nodes, Parent DFT calculations (RemoteData)
- **Output data types**: ArrayData (quasiparticle energies), BandsData

## Interfaces & Ecosystem
- **Yambo**: The backend engine
- **AiiDA-QuantumESPRESSO**: Typically used for the precursor DFT calculation
- **Materials Cloud**: Used to share GW datasets

## Workflow and Usage
1. Run DFT (QE `pw.x`).
2. Run `p2y` and `yambo -i` (initialization).
3. Submit `YamboConvergence` WorkChain to find optimal parameters.
4. Run production GW calculation.

## Performance Characteristics
- Essential for managing the complexity of GW convergence
- Handles multi-step dependencies (DFT -> p2y -> init -> run)

## Application Areas
- High-throughput GW band gaps
- Optical spectra of 2D materials
- Benchmarking MBPT methods

## Community and Support
- Developed by Yambo team and CNR-NANO (Modena)
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/yambo-code/aiida-yambo
2. Documentation: https://aiida-yambo.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: GW/BSE workflows, convergence automation
