# OpenMM

## Official Resources
- Homepage: https://openmm.org/
- Documentation: http://docs.openmm.org/
- Source Repository: https://github.com/openmm/openmm
- License: MIT/LGPL

## Overview
OpenMM is a high-performance toolkit for molecular simulation. It can be used as a library, or as an application. It provides a combination of extreme flexibility (through custom forces and integrators), openness, and high performance (especially on recent GPUs) that make it unique among simulation codes.

**Scientific domain**: Molecular dynamics, GPU-accelerated simulations, custom force fields  
**Target user community**: Researchers needing flexible, high-performance MD with Python API

## Theoretical Methods
- Classical Newtonian dynamics
- Langevin dynamics
- Brownian dynamics
- Custom integrators
- Multiple force field support (AMBER, CHARMM, OPLS)
- Custom forces via Python API
- Alchemical free energy methods

## Capabilities (CRITICAL)
- GPU-accelerated MD (CUDA, OpenCL, CPU)
- Python API for full control
- Custom forces and integrators
- Implicit and explicit solvent
- Replica exchange
- Alchemical transformations
- Coarse-grained simulations
- Extensible plugin architecture

## Key Strengths

### Flexibility:
- Custom forces via Python
- Custom integrators
- Plugin architecture
- Full API access

### Performance:
- Excellent GPU acceleration
- CUDA and OpenCL support
- Competitive with specialized codes
- Mixed precision support

## Inputs & Outputs
- **Input formats**:
  - PDB structures
  - Amber prmtop/inpcrd
  - CHARMM PSF
  - GROMACS top/gro
  - OpenMM XML
  
- **Output data types**:
  - Trajectories (DCD, PDB, XTC)
  - State files
  - Checkpoint files

## Interfaces & Ecosystem
- **Python**: Native API
- **MDTraj**: Trajectory analysis
- **OpenMMTools**: Enhanced sampling
- **YANK**: Alchemical free energy
- **OpenFF**: Open Force Field

## Advanced Features
- **Custom forces**: Define any mathematical expression
- **Custom integrators**: Implement novel algorithms
- **Alchemical methods**: Free energy perturbation
- **Drude oscillators**: Polarizable force fields
- **Constant pH**: pH-dependent simulations
- **REST2**: Replica exchange with solute tempering

## Performance Characteristics
- Excellent GPU performance
- Near-linear scaling on single GPU
- Multi-GPU support
- Mixed precision for speed

## Computational Cost
- GPU provides 100x+ speedup over CPU
- Competitive with GROMACS/AMBER on GPU
- Custom forces may reduce performance
- Overall: Excellent for GPU systems

## Best Practices
- Use CUDA platform when available
- Enable mixed precision for speed
- Use PME for electrostatics
- Validate custom forces carefully

## Limitations & Known Constraints
- Single-node focus (limited multi-node)
- Custom forces slower than built-in
- Learning curve for advanced features
- Some features GPU-only

## Application Areas
- Drug discovery
- Protein dynamics
- Method development
- Custom force field research
- Free energy calculations
- Enhanced sampling development

## Comparison with Other Codes
- **vs GROMACS**: OpenMM more flexible/customizable, GROMACS faster for standard simulations
- **vs AMBER**: OpenMM open-source with Python API, AMBER commercial with extensive force fields
- **vs NAMD**: OpenMM better GPU single-node, NAMD better multi-node scaling
- **vs LAMMPS**: OpenMM biomolecular focus, LAMMPS materials science focus
- **Unique strength**: Custom forces/integrators via Python, extreme flexibility, excellent GPU performance

## Community and Support
- Active development (Stanford)
- Large user community
- GitHub issues
- Mailing list
- Extensive documentation

## Verification & Sources
**Primary sources**:
1. Website: https://openmm.org/
2. P. Eastman et al., PLoS Comput. Biol. 13, e1005659 (2017)
3. P. Eastman et al., J. Chem. Theory Comput. 9, 461 (2013)

**Secondary sources**:
1. OpenMM tutorials and cookbooks
2. OpenMMTools documentation
3. Extensive published applications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT/LGPL)
- Academic citations: >3000
- Active development: Regular releases, 10,000+ commits
- Community: Large user base, active forums
