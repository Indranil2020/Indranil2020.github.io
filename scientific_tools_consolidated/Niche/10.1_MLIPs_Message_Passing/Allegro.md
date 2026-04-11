# Allegro

## Official Resources
- Homepage: https://github.com/mir-group/allegro
- Documentation: https://github.com/mir-group/allegro
- Source Repository: https://github.com/mir-group/allegro
- License: MIT License

## Overview
Allegro is a strictly local equivariant deep learning interatomic potential. It is built on the same principles as NequIP (E(3)-equivariance) but is designed to be strictly local (no message passing beyond a cutoff) and massively parallel. This allows it to scale to extremely large systems (millions of atoms) while maintaining high accuracy.

**Scientific domain**: Machine learning potentials, large-scale MD  
**Target user community**: Researchers simulating very large systems (proteins, cracks, grain boundaries)

## Capabilities (CRITICAL)
- **Strict Locality**: Interactions are strictly limited to a cutoff radius, enabling efficient parallelization.
- **Equivariance**: Uses tensor products for high accuracy.
- **Scalability**: Linear scaling with number of atoms, excellent strong scaling on GPUs.
- **LAMMPS**: Integration for large-scale MD.

**Sources**: Allegro GitHub, Nat. Commun. 14, 2038 (2023)

## Inputs & Outputs
- **Input formats**: Training data (extxyz)
- **Output data types**: PyTorch models

## Interfaces & Ecosystem
- **NequIP**: Share the same codebase/infrastructure.
- **LAMMPS**: Primary deployment target.

## Workflow and Usage
1. Train model using `nequip-train` with Allegro config.
2. Deploy for LAMMPS.
3. Run large-scale MD.

## Performance Characteristics
- Faster than message-passing networks for large systems.
- High parallel efficiency.

## Application Areas
- Large-scale fracture simulations
- Biomacromolecules
- Electrolytes

## Community and Support
- Developed by Kozinsky Group (Harvard)

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mir-group/allegro
2. Publication: A. Musaelian et al., Nat. Commun. 14, 2038 (2023)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Large-scale equivariant MD
