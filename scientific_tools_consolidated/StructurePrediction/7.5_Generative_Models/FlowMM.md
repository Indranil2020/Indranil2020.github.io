# FlowMM

## Overview
FlowMM is a generative model for materials using Riemannian flow matching. It generates crystal structures by learning flows on the Riemannian manifold of periodic structures, achieving state-of-the-art performance on CSP tasks.

## Theoretical Basis
- Riemannian flow matching
- Continuous normalizing flows
- Periodic structure manifold
- Equivariant architecture
- Conditional flow matching

## Key Capabilities
- Crystal structure prediction
- De novo material generation
- Riemannian geometry handling
- State-of-the-art performance
- FlowLLM extension (LLM base)

**Sources**: arXiv:2406.04713, Facebook Research

## Key Strengths

### Methodology:
- Riemannian flow matching
- Proper geometry handling
- Continuous flows

### Performance:
- State-of-the-art results
- Multiple benchmarks
- Efficient sampling

### Extensions:
- FlowLLM (LLM integration)
- CrystalLLM base
- Flexible architecture

## Inputs & Outputs
- **Input formats**: Chemical composition (CSP), nothing (de novo)
- **Output data types**: Generated crystal structures

## Interfaces & Ecosystem
- **Framework**: PyTorch, PyTorch Lightning
- **Datasets**: Standard CSP benchmarks
- **Extensions**: FlowLLM, CrystalLLM

## Workflow and Usage
1. Train on crystal dataset
2. Define generation task (CSP or de novo)
3. Sample from flow model
4. Generate structures
5. Validate candidates

## Performance Characteristics
- GPU-accelerated
- Efficient sampling
- Good generation quality

## Computational Cost
- Training: GPU hours
- Sampling: moderate
- Validation: DFT

## Best Practices
- Use appropriate training data
- Generate multiple samples
- Validate with DFT
- Consider FlowLLM for enhanced performance

## Limitations & Known Constraints
- Training data dependent
- Complex implementation
- Requires GPU

## Application Areas
- Crystal structure prediction
- De novo material generation
- Materials discovery
- Generative materials science

## Comparison with Other Codes
- **vs DiffCSP**: FlowMM flow-based, DiffCSP diffusion
- **vs CDVAE**: Different generative approach
- **Unique strength**: Riemannian flow matching, Facebook Research

## Community and Support
- Open-source (GitHub)
- Facebook Research
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/facebookresearch/flowmm
2. arXiv: https://arxiv.org/abs/2406.04713

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: ACTIVE (Meta)
- Applications: Crystal structure prediction, generative modeling
