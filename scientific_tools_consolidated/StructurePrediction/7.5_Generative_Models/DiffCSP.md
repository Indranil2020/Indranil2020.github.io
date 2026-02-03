# DiffCSP

## Overview
DiffCSP is a crystal structure prediction method using joint equivariant diffusion on lattice and atomic coordinates. It treats CSP as a conditional generation task given chemical composition.

## Theoretical Basis
- Equivariant diffusion models
- Joint lattice and coordinate diffusion
- Periodic E(3) equivariance
- Conditional generation
- Score-based denoising

## Key Capabilities
- Crystal structure prediction from composition
- Joint lattice and atom diffusion
- Equivariant architecture
- Conditional generation
- State-of-the-art performance

**Sources**: NeurIPS 2023, arXiv:2309.04475

## Key Strengths

### Equivariance:
- Periodic E(3) equivariance
- Proper symmetry handling
- Rotation/translation invariance

### Joint Diffusion:
- Lattice and coordinates together
- Consistent generation
- Physical constraints

### Performance:
- State-of-the-art results
- Multiple benchmarks
- Competitive accuracy

## Inputs & Outputs
- **Input formats**: Chemical composition, atom types
- **Output data types**: Predicted crystal structures

## Interfaces & Ecosystem
- **Framework**: PyTorch
- **Datasets**: MP-20, Perov-5, Carbon-24, MPTS-52
- **Evaluation**: Match rate, RMSD

## Workflow and Usage
1. Train on crystal structure dataset
2. Input chemical composition
3. Run diffusion sampling
4. Generate candidate structures
5. Rank and validate

## Performance Characteristics
- GPU-accelerated
- Fast sampling
- Good match rates

## Computational Cost
- Training: GPU hours
- Sampling: moderate
- Validation: DFT-dependent

## Best Practices
- Use appropriate training data
- Generate multiple samples
- Validate top candidates
- Check structural validity

## Limitations & Known Constraints
- Training data dependent
- Complex compositions challenging
- Requires validation

## Application Areas
- Crystal structure prediction
- Materials discovery
- Composition-to-structure
- High-throughput screening

## Comparison with Other Codes
- **vs CDVAE**: DiffCSP pure diffusion, CDVAE VAE+diffusion
- **vs EquiCSP**: Similar approach, different implementation
- **Unique strength**: Joint equivariant diffusion, NeurIPS 2023

## Community and Support
- Open-source (GitHub)
- Academic development
- Active research area

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/jiaor17/DiffCSP
2. Paper: NeurIPS 2023
3. arXiv: https://arxiv.org/abs/2309.04475

**Secondary sources**:
1. DiffCSP++ (ICLR 2024): https://github.com/jiaor17/DiffCSP-PP

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: ACTIVE
- Applications: Crystal structure prediction, generative modeling
