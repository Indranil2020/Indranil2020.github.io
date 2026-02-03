# CDVAE

## Overview
CDVAE (Crystal Diffusion Variational Autoencoder) is an SE(3)-invariant autoencoder for generating periodic crystal structures. It combines variational autoencoders with diffusion models to generate stable materials.

## Theoretical Basis
- Variational autoencoder (VAE)
- Score-based diffusion model
- SE(3) equivariance
- Periodic structure handling
- Energy-guided generation

## Key Capabilities
- Periodic material generation
- Stable structure generation
- Property optimization
- Reconstruction and generation
- SE(3) invariance

**Sources**: ICLR 2022, arXiv:2110.06197

## Key Strengths

### Architecture:
- SE(3)-invariant
- Diffusion-based decoder
- VAE latent space

### Generation:
- Stable structures
- Diverse compositions
- Property-guided

### Validation:
- Energy evaluation
- Stability checks
- Benchmark datasets

## Inputs & Outputs
- **Input formats**: Crystal structures (training), latent vectors (generation)
- **Output data types**: Generated crystal structures

## Interfaces & Ecosystem
- **Framework**: PyTorch, PyTorch Geometric
- **Datasets**: MP-20, Carbon-24, Perov-5
- **Evaluation**: DFT validation

## Workflow and Usage
1. Train on crystal structure dataset
2. Encode structures to latent space
3. Sample from latent space
4. Decode with diffusion process
5. Evaluate generated structures

## Performance Characteristics
- GPU-accelerated training
- Fast generation after training
- Good reconstruction accuracy

## Computational Cost
- Training: GPU hours
- Generation: fast
- Evaluation: DFT-dependent

## Best Practices
- Use appropriate training data
- Validate with DFT
- Check stability metrics
- Use property guidance

## Limitations & Known Constraints
- Training data dependent
- May generate unstable structures
- Limited to training distribution

## Application Areas
- Materials discovery
- Crystal structure generation
- Property-guided design
- Stable material prediction

## Comparison with Other Codes
- **vs DiffCSP**: CDVAE VAE-based, DiffCSP pure diffusion
- **vs MatterGen**: Different architectures
- **Unique strength**: VAE+diffusion, SE(3) invariance, ICLR 2022 landmark

## Community and Support
- Open-source (GitHub)
- MIT/Meta AI development
- Widely cited

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/txie-93/cdvae
2. Paper: ICLR 2022
3. arXiv: https://arxiv.org/abs/2110.06197

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: LANDMARK (2022)
- Applications: Crystal structure generation, materials discovery
