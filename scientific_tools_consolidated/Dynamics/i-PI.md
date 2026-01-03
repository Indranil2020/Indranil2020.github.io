# i-PI

## Official Resources
- Homepage: http://ipi-code.org/
- Documentation: http://ipi-code.org/documentation/
- Source Repository: https://github.com/i-pi/i-pi
- License: GNU General Public License v3.0

## Overview
i-PI is a universal force engine interface that decouples the evolution of nuclear coordinates from the evaluation of the potential energy surface. It acts as a client-server driver for molecular dynamics, enabling advanced quantum nuclear effects (path integral MD), thermodynamic integration, and enhanced sampling techniques using any electronic structure code that supports the i-PI socket interface.

**Scientific domain**: Ab-initio molecular dynamics, path integral MD, quantum nuclear effects  
**Target user community**: Computational chemists, materials scientists, quantum dynamics researchers

## Theoretical Methods
- Path Integral Molecular Dynamics (PIMD)
- Ring Polymer Molecular Dynamics (RPMD)
- Centroid Molecular Dynamics (CMD)
- Thermostatted Ring Polymer MD (TRPMD)
- Replica exchange MD
- Thermodynamic integration
- Metadynamics (via PLUMED)
- Generalized Langevin Equation (GLE) thermostats

## Capabilities (CRITICAL)
- Decoupled force evaluation (client-server model)
- Quantum nuclear effects via path integrals
- Advanced thermostats (GLE, pile)
- Multiple time step integration
- Geometry optimization
- Phonon calculations
- Isotope fractionation
- Integration with virtually any DFT/force code (VASP, QE, CP2K, LAMMPS, etc.)
- Python-based driver with C++ clients

**Sources**: i-PI documentation, Comp. Phys. Comm. 205, 106 (2016)

## Inputs & Outputs
- **Input formats**: XML input file controlling dynamics, sockets, ensembles
- **Output data types**: Trajectories (pdb, xyz), properties, restart files

## Interfaces & Ecosystem
- **Clients (Force Engines)**: VASP, Quantum ESPRESSO, CP2K, LAMMPS, FHI-aims, Siesta, DFTB+, xTB, and many others
- **PLUMED**: Interface for enhanced sampling
- **ASE**: Compatible via calculators

## Workflow and Usage
1. Start i-PI server: `i-pi input.xml`
2. Start force engine client(s): `lmp_mpi -in in.lammps` (configured for i-PI)
3. i-PI sends positions to client
4. Client calculates forces and energy, returns to i-PI
5. i-PI propagates dynamics

## Performance Characteristics
- Minimal overhead from Python driver
- Parallelism via multiple force clients
- Efficient socket communication

## Application Areas
- Water and aqueous solutions (nuclear effects)
- Hydrogen storage materials
- Proton transfer reactions
- Low-temperature dynamics
- Isotope effects in materials

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Active development (Michele Ceriotti group)
- User forum/mailing list

## Verification & Sources
**Primary sources**:
1. Homepage: http://ipi-code.org/
2. GitHub: https://github.com/i-pi/i-pi
3. Publication: Comp. Phys. Comm. 205, 106 (2016)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (EPFL)
- Applications: Universal force engine, path integrals, quantum nuclear effects
