# AmpTorch

## Official Resources
- Homepage: https://github.com/ulissigroup/amptorch
- Documentation: https://amptorch.readthedocs.io/
- Source Repository: https://github.com/ulissigroup/amptorch
- License: Apache-2.0

## Overview
AmpTorch is a PyTorch implementation of the Atomistic Machine-learning Package (AMP) for training neural network potentials. It provides GPU acceleration and modern deep learning features for developing interatomic potentials.

**Scientific domain**: Neural network potentials, GPU-accelerated training  
**Target user community**: Researchers training custom ML potentials

## Theoretical Methods
- Behler-Parrinello symmetry functions
- Neural network potentials
- Gaussian descriptor functions
- PyTorch automatic differentiation

## Capabilities (CRITICAL)
- GPU-accelerated training
- PyTorch backend
- Multiple descriptor types
- ASE calculator interface
- Custom architectures
- Transfer learning

## Key Strengths

### PyTorch Integration:
- GPU acceleration
- Modern deep learning
- Automatic differentiation
- Easy customization

### Flexibility:
- Custom descriptors
- Custom architectures
- Transfer learning

## Inputs & Outputs
- **Input formats**:
  - ASE trajectory files
  - VASP OUTCAR
  
- **Output data types**:
  - Energies
  - Forces
  - Model checkpoints

## Interfaces & Ecosystem
- **ASE**: Calculator
- **PyTorch**: Backend
- **AMP**: Original framework

## Advanced Features
- **GPU training**: CUDA acceleration
- **Custom descriptors**: Flexible features
- **Transfer learning**: Fine-tuning
- **Ensemble**: Multiple models

## Performance Characteristics
- Fast GPU training
- Efficient inference
- Good scaling

## Computational Cost
- Training: Hours (GPU)
- Inference: Fast
- Overall: Efficient

## Best Practices
- Use GPU for training
- Validate on test set
- Use appropriate descriptors

## Limitations & Known Constraints
- Requires training data
- Descriptor choice important
- Less active than alternatives

## Application Areas
- Catalysis
- Surface science
- Materials modeling

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ulissigroup/amptorch

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, Apache-2.0)
