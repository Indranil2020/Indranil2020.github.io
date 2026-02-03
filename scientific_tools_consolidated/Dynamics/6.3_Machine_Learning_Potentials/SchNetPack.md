# SchNetPack

## Official Resources
- Homepage: https://schnetpack.readthedocs.io/
- Documentation: https://schnetpack.readthedocs.io/
- Source Repository: https://github.com/atomistic-machine-learning/schnetpack
- License: MIT

## Overview
SchNetPack is a toolbox for deep neural networks in atomistic systems. It provides implementations of SchNet, PaiNN, and other equivariant neural network architectures for learning potential energy surfaces and molecular properties.

**Scientific domain**: Deep learning for molecules, neural network potentials, property prediction  
**Target user community**: Researchers developing and applying neural network potentials

## Theoretical Methods
- SchNet continuous-filter convolutions
- PaiNN equivariant message passing
- E(3)-equivariant architectures
- Property prediction
- Potential energy surfaces

## Capabilities (CRITICAL)
- SchNet architecture
- PaiNN equivariant networks
- Property prediction
- Force field training
- ASE calculator interface
- PyTorch Lightning training

## Key Strengths

### Architectures:
- SchNet (original)
- PaiNN (equivariant)
- Modular design
- Easy customization

### Training:
- PyTorch Lightning
- Multi-GPU support
- Extensive logging

## Inputs & Outputs
- **Input formats**: ASE databases, XYZ, QM datasets
- **Output data types**: Energies, forces, properties, models

## Interfaces & Ecosystem
- **ASE**: Calculator interface
- **PyTorch**: Backend
- **PyTorch Lightning**: Training

## Advanced Features
- **SchNet**: Continuous-filter convolutions
- **PaiNN**: Polarizable equivariant
- **Property prediction**: Multiple outputs
- **Transfer learning**: Fine-tuning support

## Performance Characteristics
- Good training speed
- GPU acceleration
- Efficient inference

## Computational Cost
- Training: Hours to days (GPU)
- Inference: Fast
- Overall: Efficient

## Best Practices
- Use PaiNN for accuracy
- Validate on test set
- Use appropriate cutoffs

## Limitations & Known Constraints
- Requires training data
- Architecture choice matters
- GPU recommended

## Application Areas
- Molecular property prediction
- Potential energy surfaces
- Drug discovery
- Materials science

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/atomistic-machine-learning/schnetpack
2. K.T. Schütt et al., J. Chem. Theory Comput. 15, 448 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
