# AMBER (Assisted Model Building with Energy Refinement)

## Official Resources
- Homepage: https://ambermd.org/
- Documentation: https://ambermd.org/doc12/
- Source Repository: https://gitlab.com/amber-md/amber (AmberTools is open)
- License: Mixed (AmberTools: GPL/LGPL, AMBER: Commercial/Academic license)

## Overview
AMBER refers to two things: a set of molecular mechanical force fields for the simulation of biomolecules (which are in the public domain), and a package of molecular simulation programs. The software package includes AmberTools (open source) for setup and analysis, and the AMBER MD engine (pmemd) which is highly optimized for GPU acceleration.

**Scientific domain**: Biomolecular simulation, drug design, GPU molecular dynamics  
**Target user community**: Computational chemists, structural biologists, pharmaceutical researchers

## Theoretical Methods
- Classical Molecular Dynamics (PMEMD)
- Generalized Born implicit solvent models
- PME (Particle Mesh Ewald) for electrostatics
- TI (Thermodynamic Integration) for free energy
- QM/MM (Quantum Mechanics/Molecular Mechanics)
- pH-constant molecular dynamics
- Gaussian Accelerated Molecular Dynamics (GaMD)
- Nudged Elastic Band (NEB)

## Capabilities (CRITICAL)
- GPU-accelerated MD (among the fastest available)
- Simulation of proteins, nucleic acids, carbohydrates
- Comprehensive analysis tools (cpptraj)
- NMR refinement structure calculation
- Implicit solvent simulations (GB/PB)
- Force field development (GAFF, ff14SB, OL15, etc.)
- Free energy calculations (TI, MM-PBSA/GBSA)

**Sources**: AMBER website, J. Chem. Inf. Model. 58, 2043 (2018)

## Key Strengths

### Force Fields:
- AMBER force fields (ff14SB, ff19SB)
- GAFF for small molecules
- Extensively validated
- Regular updates

### GPU Performance:
- Industry-leading pmemd.cuda
- Multi-GPU support
- Long timescales routine

### Analysis:
- cpptraj comprehensive analysis
- pytraj Python bindings
- MM-PBSA/GBSA

## Inputs & Outputs
- **Input formats**: prmtop (topology/params), inpcrd (coordinates), mdin (control parameters)
- **Output data types**: mdcrd/NetCDF (trajectory), mdout (log), rst (restart)

## Interfaces & Ecosystem
- **AmberTools**: Essential for setup (tleap, antechamber) and analysis (cpptraj)
- **VMD/Chimera**: Visualization
- **Python**: pytraj (Python bindings for cpptraj)
- **PLUMED**: Interface available

## Workflow and Usage
1. Prepare structure: `pdb4amber -i protein.pdb`
2. Parameterize ligands: `antechamber` (GAFF)
3. Build system: `tleap` (solvate, add ions, save prmtop/inpcrd)
4. Energy minimization: `pmemd -O -i min.in ...`
5. Heating/Equilibration: `pmemd.cuda -O -i heat.in ...`
6. Production: `pmemd.cuda -O -i prod.in ...`
7. Analysis: `cpptraj`

## Performance Characteristics
- **GPU**: Industry-leading GPU performance (pmemd.cuda)
- **Scaling**: Excellent on single/multi-GPU nodes
- **Efficiency**: Optimized for long timescales

## Computational Cost
- Excellent GPU performance
- Microseconds/day on modern GPUs
- Efficient for biomolecules
- Overall: Industry-leading for drug discovery

## Best Practices
- Use tleap for system setup
- Validate force field choice
- Use cpptraj for analysis
- Enable GPU acceleration
- Check for clashes before production

## Limitations & Known Constraints
- Commercial license for pmemd
- Less flexible than LAMMPS
- Biomolecular focus
- Complex parameter files

## Application Areas
- Drug discovery (lead optimization)
- Protein-ligand binding
- Conformational sampling
- Nucleic acid dynamics
- Refinement of NMR/X-ray structures

## Comparison with Other Codes
- **vs GROMACS**: AMBER better GPU single-node, GROMACS open-source
- **vs NAMD**: AMBER faster GPU, NAMD better multi-node
- **vs CHARMM**: AMBER more GPU-focused, CHARMM more methods
- **Unique strength**: AMBER force fields, GPU performance, drug discovery focus

## Community and Support
- AmberTools: Open source community
- AMBER: Licensed software with support
- Active mailing list (amber@ambermd.org)
- Annual workshops

## Verification & Sources
**Primary sources**:
1. Homepage: https://ambermd.org/
2. AmberTools: https://ambermd.org/AmberTools.php
3. Publication: Case et al., J. Comput. Chem. 26, 1668 (2005)

**Secondary sources**:
1. AMBER tutorials
2. cpptraj documentation
3. Extensive published applications (>30,000 citations)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: MIXED (AmberTools open, PMEMD licensed)
- Development: ACTIVE (Rutgers, UCSF, etc.)
- Applications: Force fields, GPU MD, drug design, biomolecules
