# NAMD (Nanoscale Molecular Dynamics)

## Official Resources
- Homepage: https://www.ks.uiuc.edu/Research/namd/
- Documentation: https://www.ks.uiuc.edu/Research/namd/current/ug/
- Source Repository: https://gitlab.com/tcbg/namd (and UIUC download)
- License: UIUC Open Source License (Non-exclusive, simplified)

## Overview
NAMD is a parallel molecular dynamics code designed for high-performance simulation of large biomolecular systems. Developed by the Theoretical and Computational Biophysics Group (TCBG) at UIUC, it is renowned for its scalability on massive supercomputers and its ease of use. NAMD is file-compatible with AMBER, CHARMM, and X-PLOR.

**Scientific domain**: Biomolecular simulation, massive parallel MD
**Target user community**: Biophysicists, biochemists, supercomputing users

## Theoretical Methods
- Classical Molecular Dynamics
- Steered Molecular Dynamics (SMD)
- Interactive Molecular Dynamics (IMD)
- Alchemical Free Energy Perturbation (FEP)
- Replica Exchange MD
- Constant pH MD
- QM/MM Simulation
- Multiple Time Stepping

## Capabilities (CRITICAL)
- Extreme scalability (up to hundreds of thousands of cores)
- GPU acceleration (CUDA)
- Tcl-based scripting interface
- Compatibility with CHARMM, AMBER, and PLUMED
- Interface with VMD for visualization and setup
- QwikMD for easy setup
- Grid-based forces (GridForces) for multiscale modeling

**Sources**: NAMD website, J. Comput. Chem. 26, 1781 (2005)

## Inputs & Outputs
- **Input formats**: Configuration file (.conf), PDB/PSF files (structure/topology), Parameter files
- **Output data types**: DCD trajectories, Output logs, Restart files, Velocity files

## Interfaces & Ecosystem
- **VMD**: Deeply integrated visualization and analysis tool
- **QwikMD**: Easy-to-use GUI within VMD
- **PLUMED**: Supported via patching
- **Colvars**: Collective variables module (native)

## Workflow and Usage
1. **Setup**: Use VMD (PSFgen) or QwikMD to build the system (PSF and PDB).
2. **Configuration**: Create a NAMD config file specifying parameters, force fields, and run steps.
3. **Run**: `namd2 +p<cores> config.conf > output.log`
4. **Analysis**: Load DCD and PSF into VMD.

## Performance Characteristics
- **Scalability**: One of the most scalable MD codes available (Charm++ runtime)
- **GPU**: Excellent CUDA performance
- **Efficiency**: Optimized for large systems (millions of atoms)

## Application Areas
- Viral capsids and large assemblies
- Membrane proteins and channels
- Protein folding
- Free energy calculations
- Mechanobiology (Steered MD)

## Community and Support
- Developed by TCBG at UIUC (NIH funded)
- Extensive tutorials and wiki
- Very active mailing list
- Free for academic use

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.ks.uiuc.edu/Research/namd/
2. Publication: Phillips et al., J. Comput. Chem. 26, 1781 (2005)
3. Repository: https://gitlab.com/tcbg/namd

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN
- Development: ACTIVE (UIUC)
- Applications: High-performance biomolecular MD
