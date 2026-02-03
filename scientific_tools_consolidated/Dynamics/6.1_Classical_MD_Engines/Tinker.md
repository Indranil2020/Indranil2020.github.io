# Tinker

## Official Resources
- Homepage: https://dasher.wustl.edu/tinker/
- Documentation: https://dasher.wustl.edu/tinker/distribution/doc/
- Source Repository: https://github.com/TinkerTools/tinker
- License: Custom (free for non-commercial)

## Overview
Tinker is a complete and general package for molecular mechanics and dynamics, with special features for biopolymers. It is particularly known for its implementation of polarizable force fields, especially the AMOEBA (Atomic Multipole Optimized Energetics for Biomolecular Applications) force field.

**Scientific domain**: Polarizable force fields, biomolecular simulations, AMOEBA  
**Target user community**: Researchers using polarizable force fields

## Theoretical Methods
- Classical molecular dynamics
- Polarizable force fields (AMOEBA)
- Atomic multipoles
- Induced dipoles
- Multiple force fields (AMBER, CHARMM, OPLS, MM3)
- Free energy perturbation

## Capabilities (CRITICAL)
- AMOEBA polarizable force field
- Atomic multipole electrostatics
- Multiple force field support
- Free energy calculations
- Normal mode analysis
- Geometry optimization
- GPU acceleration (Tinker-HP)

## Key Strengths

### Polarizable Force Fields:
- AMOEBA implementation
- Atomic multipoles
- Induced dipoles
- Accurate electrostatics

### Versatility:
- Many force fields
- Analysis tools
- Optimization methods

## Inputs & Outputs
- **Input formats**:
  - XYZ coordinates
  - PDB structures
  - Tinker key files
  
- **Output data types**:
  - Trajectories
  - Energy files
  - Analysis output

## Interfaces & Ecosystem
- **Tinker-HP**: GPU/parallel version
- **Tinker-OpenMM**: OpenMM interface
- **Poltype**: Parameterization tool
- **FFX**: Force Field X

## Advanced Features
- **AMOEBA**: Advanced polarizable FF
- **Multipoles**: Higher-order electrostatics
- **Tinker-HP**: Massively parallel version
- **QM/MM**: Quantum mechanics coupling
- **Free energy**: FEP and TI methods
- **Optimization**: Various minimizers

## Performance Characteristics
- Tinker-HP: Excellent parallel scaling
- GPU support via Tinker-HP
- Efficient for polarizable FF
- Good for medium systems

## Computational Cost
- Polarizable FF more expensive than fixed-charge
- Tinker-HP provides speedup
- GPU acceleration available
- Overall: Moderate (polarizable overhead)

## Best Practices
- Use Tinker-HP for large systems
- Validate AMOEBA parameters
- Use appropriate polarization damping
- Check energy conservation

## Limitations & Known Constraints
- Polarizable FF more expensive
- Steeper learning curve
- Less community than GROMACS/AMBER
- Parameter development needed

## Application Areas
- Polarizable simulations
- Protein-ligand binding
- Ion solvation
- Accurate electrostatics
- Force field development

## Community and Support
- Active development
- Mailing list
- Documentation
- Tutorials available

## Verification & Sources
**Primary sources**:
1. Website: https://dasher.wustl.edu/tinker/
2. J.W. Ponder et al., J. Phys. Chem. B 114, 2549 (2010)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: Available (GitHub)
- Well-documented
