# OpenMD

## Official Resources
- Homepage: https://openmd.org/
- Documentation: https://openmd.org/documentation/
- Source Repository: https://github.com/OpenMD/OpenMD
- License: BSD 3-Clause License

## Overview
OpenMD is an open source molecular dynamics engine written in C++ that is designed to simulate liquids, proteins, nanoparticles, interfaces, and other complex systems. It focuses on versatility and ease of use, with a particular emphasis on handling non-standard potentials and rigid body dynamics.

**Scientific domain**: Molecular dynamics, soft matter, metallic nanoparticles, interfaces  
**Target user community**: Physical chemists, materials scientists

## Theoretical Methods
- Classical Molecular Dynamics (NVE, NVT, NPT)
- Rigid Body Dynamics (quaternions)
- Electrostatics (SPME, damped shifted force)
- Minimization (SD, CG)
- Z-constraint methods (for free energy profiles)
- Thermodynamic Integration
- Langevin Hull method (NPT)

## Capabilities (CRITICAL)
- Simulation of atomistic and rigid-body systems
- Embedded Atom Method (EAM) for metals
- Transition metal oxides and water
- Nanoparticle simulations (melting, interfaces)
- Slab geometry electrostatics
- Fluctuating charge models (EAM-µ, electronegativity equalization)
- Restraints and external fields
- Parallel execution (MPI)

**Sources**: OpenMD website, J. Chem. Phys. 124, 024109 (2006)

## Inputs & Outputs
- **Input formats**: .omd file (XML-like structure + coordinates)
- **Output data types**: .dump (trajectory), .stat (thermodynamics), .eor (end of run)

## Interfaces & Ecosystem
- **Python**: Analysis scripts
- **VMD**: Visualization support
- **Standalone**: All-in-one input format

## Workflow and Usage
1. Create `.omd` file: Contains force field, topology, and initial coordinates
2. Run simulation: `openmd system.omd`
3. Convert output: `Dump2XYZ` or `Dump2PDB`
4. Analysis: `StaticProps`, `DynamicProps` tools

## Performance Characteristics
- Good scaling on moderate clusters
- Specialized for rigid bodies and metals
- Not as ultra-optimized as GROMACS for biomolecules but efficient for general chemistry

## Application Areas
- Metallic nanoparticles and alloys
- Lipid bilayers
- Water interfaces
- Zeolites and minerals
- Ionic liquids

## Community and Support
- Open-source (BSD)
- Developed at University of Notre Dame (Gezelter group)
- Mailing list and issue tracker

## Verification & Sources
**Primary sources**:
1. Homepage: https://openmd.org/
2. GitHub: https://github.com/OpenMD/OpenMD
3. Publication: M. A. Meineke et al., J. Comput. Chem. 26, 252 (2005)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (Notre Dame)
- Applications: MD, rigid bodies, metals, fluctuating charge
