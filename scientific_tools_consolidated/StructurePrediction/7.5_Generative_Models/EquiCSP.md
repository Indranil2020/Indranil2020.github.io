# EquiCSP

## Overview
EquiCSP (Equivariant Diffusion for Crystal Structure Prediction) is a symmetry-aware deep learning model that ensures permutation, rotation, and periodic translation equivariance during the diffusion process for crystal structure prediction.

## Theoretical Basis
- Equivariant diffusion models
- Permutation equivariance
- Rotation equivariance
- Periodic translation equivariance
- Conditional generation

## Key Capabilities
- Equivariant crystal structure prediction
- Symmetry-aware generation
- Conditional on composition
- ICML 2024 publication
- State-of-the-art performance

**Sources**: ICML 2024, arXiv:2512.07289

## Key Strengths

### Equivariance:
- Full symmetry handling
- Permutation invariance
- Rotation/translation equivariance

### Architecture:
- Diffusion-based
- Symmetry-aware design
- Proper periodic handling

### Performance:
- Competitive results
- ICML publication
- Well-validated

## Inputs & Outputs
- **Input formats**: Chemical composition
- **Output data types**: Predicted crystal structures

## Interfaces & Ecosystem
- **Framework**: PyTorch
- **Base**: Built on CDVAE/DiffCSP
- **Datasets**: Standard benchmarks

## Workflow and Usage
1. Train on crystal dataset
2. Input composition
3. Run equivariant diffusion
4. Generate structures
5. Validate candidates

## Performance Characteristics
- GPU-accelerated
- Efficient sampling
- Good accuracy

## Computational Cost
- Training: GPU hours
- Sampling: moderate
- Validation: DFT

## Best Practices
- Use appropriate training data
- Generate multiple samples
- Validate with DFT
- Check symmetry preservation

## Limitations & Known Constraints
- Training data dependent
- Complex implementation
- Requires GPU

## Application Areas
- Crystal structure prediction
- Symmetry-aware generation
- Materials discovery
- Composition-to-structure

## Comparison with Other Codes
- **vs DiffCSP**: EquiCSP more equivariant
- **vs CDVAE**: Different equivariance approach
- **Unique strength**: Full equivariance, ICML 2024

## Community and Support
- Open-source (GitHub)
- Academic development
- ICML publication

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/EmperorJia/EquiCSP
2. Paper: ICML 2024
3. arXiv: https://arxiv.org/abs/2512.07289

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: ACTIVE
- Applications: Equivariant CSP
