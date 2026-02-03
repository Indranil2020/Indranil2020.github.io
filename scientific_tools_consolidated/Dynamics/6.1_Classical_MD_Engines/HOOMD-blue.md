# HOOMD-blue

## Official Resources
- Homepage: https://glotzerlab.engin.umich.edu/hoomd-blue/
- Documentation: https://hoomd-blue.readthedocs.io/
- Source Repository: https://github.com/glotzerlab/hoomd-blue
- License: BSD-3-Clause

## Overview
HOOMD-blue is a general-purpose particle simulation toolkit optimized for performance on both GPUs and CPUs. It is designed for soft matter research but is general enough for many types of particle simulations including molecular dynamics and Monte Carlo.

**Scientific domain**: Soft matter physics, colloidal systems, polymers, Monte Carlo  
**Target user community**: Soft matter researchers, materials scientists

## Theoretical Methods
- Classical molecular dynamics
- Monte Carlo methods (NVT, NPT, Gibbs ensemble)
- Hard particle simulations
- Dissipative particle dynamics (DPD)
- Brownian dynamics
- Rigid body dynamics

## Capabilities (CRITICAL)
- GPU-accelerated MD and MC
- Hard particle Monte Carlo
- Anisotropic particles
- Rigid body dynamics
- Dissipative particle dynamics
- Custom pair potentials
- Python scripting interface

## Key Strengths

### Soft Matter Focus:
- Anisotropic particles
- Hard particle MC
- Coarse-grained models
- Polymer simulations

### GPU Performance:
- Excellent GPU scaling
- CUDA optimized
- Large system support
- Efficient neighbor lists

## Inputs & Outputs
- **Input formats**:
  - GSD files (native)
  - Python initialization
  
- **Output data types**:
  - GSD trajectories
  - Log files
  - Custom outputs

## Interfaces & Ecosystem
- **Python**: Native API
- **freud**: Analysis library
- **signac**: Workflow management
- **GSD**: File format library
- **fresnel**: Visualization

## Advanced Features
- **HPMC**: Hard particle Monte Carlo
- **DPD**: Dissipative particle dynamics
- **Rigid bodies**: Composite particles
- **Custom forces**: User-defined potentials
- **Alchemical**: Free energy methods
- **Active matter**: Self-propelled particles

## Performance Characteristics
- Excellent GPU performance
- Scales to millions of particles
- Efficient for soft matter
- Good weak scaling

## Computational Cost
- GPU provides major speedup
- Efficient for short-range forces
- MC methods very efficient
- Overall: Excellent for soft matter

## Best Practices
- Use GSD format for I/O
- Enable GPU when available
- Use appropriate neighbor list settings
- Validate with known systems

## Limitations & Known Constraints
- Less biomolecular focus
- Limited long-range electrostatics
- Specialized for soft matter
- Learning curve for MC methods

## Application Areas
- Colloidal self-assembly
- Polymer physics
- Liquid crystals
- Active matter
- Granular materials
- Nanoparticle systems

## Comparison with Other Codes
- **vs LAMMPS**: HOOMD-blue better for soft matter/MC, LAMMPS more general-purpose
- **vs ESPResSo**: HOOMD-blue stronger GPU/MC, ESPResSo better for charged systems
- **vs GROMACS**: HOOMD-blue soft matter focus, GROMACS biomolecular focus
- **Unique strength**: Hard particle Monte Carlo, anisotropic particles, excellent GPU performance, Python-native

## Community and Support
- Active development (Glotzer group, Michigan)
- GitHub issues
- Documentation
- Tutorials available

## Verification & Sources
**Primary sources**:
1. Website: https://glotzerlab.engin.umich.edu/hoomd-blue/
2. J.A. Anderson et al., Comput. Mater. Sci. 173, 109363 (2020)
3. J.A. Anderson et al., J. Comput. Phys. 227, 5342 (2008)

**Secondary sources**:
1. HOOMD-blue tutorials
2. freud analysis library documentation
3. Published soft matter applications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, BSD-3)
- Academic citations: >1500
- Active development: Regular releases
- Community: Glotzer group, soft matter community
