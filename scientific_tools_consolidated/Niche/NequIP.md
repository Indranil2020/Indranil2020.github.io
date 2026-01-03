# NequIP (Neural Equivariant Interatomic Potentials)

## Official Resources
- Homepage: https://github.com/mir-group/nequip
- Documentation: https://github.com/mir-group/nequip (README/Wiki)
- Source Repository: https://github.com/mir-group/nequip
- License: MIT License

## Overview
NequIP is a code for building E(3)-equivariant neural network interatomic potentials. It uses the `e3nn` library to ensure that the learned potentials respect rotation and translation symmetries and parity by construction. This data efficiency allows NequIP to achieve high accuracy with very small training sets compared to invariant models.

**Scientific domain**: Machine learning potentials, equivariant neural networks  
**Target user community**: MD users, ML researchers

## Capabilities (CRITICAL)
- **Equivariance**: E(3)-equivariant features (vectors, tensors) used throughout the network.
- **Data Efficiency**: High accuracy with few (100-1000) training structures.
- **LAMMPS**: Interface for running MD with NequIP potentials via `pair_nequip`.
- **ASE**: ASE calculator interface.

**Sources**: NequIP GitHub, Nat. Commun. 13, 2453 (2022)

## Inputs & Outputs
- **Input formats**: Extended XYZ (for training), YAML config
- **Output data types**: PyTorch model (.pth), deployed model for LAMMPS

## Interfaces & Ecosystem
- **PyTorch**: Core framework.
- **e3nn**: Library for equivariant operations.
- **LAMMPS**: Plugin available.
- **ASE**: Integration.

## Workflow and Usage
1. Prepare training data (extxyz with energy/forces).
2. Create config YAML.
3. Train: `nequip-train config.yaml`
4. Deploy: `nequip-deploy build --train-dir results/ model.pth`
5. Run LAMMPS MD.

## Performance Characteristics
- Slower inference than simple invariant models (due to tensor products).
- Extremely high accuracy and stability.
- Excellent for complex materials where angular dependence is critical.

## Application Areas
- Phase transitions
- Reaction dynamics
- Complex oxides
- Liquid structures

## Community and Support
- Developed by Kozinsky Group (Harvard)
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mir-group/nequip
2. Publication: S. Batzner et al., Nat. Commun. 13, 2453 (2022)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Equivariant ML potentials
