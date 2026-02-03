# SevenNet

## Official Resources
- Homepage: https://github.com/MDIL-SNU/SevenNet
- Documentation: https://sevennet.readthedocs.io/
- Source Repository: https://github.com/MDIL-SNU/SevenNet
- License: GPL-3.0

## Overview
SevenNet is a graph neural network interatomic potential package that supports efficient multi-GPU parallel molecular dynamics simulations. It provides pretrained universal potentials and enables large-scale simulations with excellent parallel scaling.

**Scientific domain**: Scalable ML potentials, parallel MD, graph neural networks  
**Target user community**: Researchers needing scalable ML potentials for large systems

## Theoretical Methods
- Graph neural networks
- E(3)-equivariant features
- Multi-GPU parallelization
- Universal potential training

## Capabilities (CRITICAL)
- Multi-GPU parallel MD
- Pretrained universal models
- LAMMPS integration
- Large-scale simulations
- Efficient scaling
- ASE calculator

## Key Strengths

### Parallel Scaling:
- Multi-GPU support
- Efficient parallelization
- Large systems
- Good weak scaling

### Pretrained Models:
- Universal potentials
- Ready to use
- Fine-tuning support

## Inputs & Outputs
- **Input formats**:
  - ASE Atoms
  - LAMMPS data
  
- **Output data types**:
  - Energies
  - Forces
  - Stresses
  - Model files

## Interfaces & Ecosystem
- **LAMMPS**: pair_style
- **ASE**: Calculator
- **PyTorch**: Backend

## Advanced Features
- **Multi-GPU**: Parallel inference
- **Universal models**: Pretrained
- **LAMMPS**: Large-scale MD
- **Fine-tuning**: Transfer learning

## Performance Characteristics
- Excellent parallel scaling
- Multi-GPU efficient
- Fast inference
- Good for large systems

## Computational Cost
- Pretrained: Ready to use
- Inference: Fast, scales well
- Training: GPU hours
- Overall: Excellent scaling

## Best Practices
- Use multi-GPU for large systems
- Start with pretrained models
- Validate carefully
- Fine-tune if needed

## Limitations & Known Constraints
- Requires multiple GPUs for best performance
- Active development
- Documentation evolving

## Application Areas
- Large-scale materials simulations
- Parallel MD
- Materials discovery
- High-throughput screening

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/MDIL-SNU/SevenNet
2. Y. Park et al., J. Chem. Theory Comput. (2024)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-3.0)
