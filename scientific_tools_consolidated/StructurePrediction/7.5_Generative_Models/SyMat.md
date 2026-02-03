# SyMat

## Overview
SyMat (Symmetry-aware generation of periodic Materials) is a generative model that incorporates crystal symmetry into the generation process. It uses a symmetry-aware probabilistic model in the coordinate diffusion process.

## Theoretical Basis
- Symmetry-aware generation
- Score-based diffusion
- Space group handling
- Wyckoff position awareness
- Periodic structure generation

## Key Capabilities
- Symmetry-aware crystal generation
- Space group preservation
- Random generation
- Property optimization
- NeurIPS 2023 publication

**Sources**: NeurIPS 2023, arXiv:2307.02707

## Key Strengths

### Symmetry:
- Space group awareness
- Wyckoff positions
- Symmetry preservation

### Generation:
- Random generation
- Property optimization
- Diverse structures

### Theory:
- Invariant to symmetry transformations
- Proper periodic handling
- Well-founded

## Inputs & Outputs
- **Input formats**: Space group (optional), composition
- **Output data types**: Generated crystal structures

## Interfaces & Ecosystem
- **Framework**: PyTorch
- **Datasets**: Standard benchmarks
- **Evaluation**: Stability metrics

## Workflow and Usage
1. Train on crystal dataset
2. Specify constraints (optional)
3. Run symmetry-aware diffusion
4. Generate structures
5. Validate candidates

## Performance Characteristics
- GPU-accelerated
- Symmetry-preserving
- Good generation quality

## Computational Cost
- Training: GPU hours
- Sampling: moderate
- Validation: DFT

## Best Practices
- Use symmetry constraints when known
- Generate diverse samples
- Validate with DFT
- Check symmetry preservation

## Limitations & Known Constraints
- Training data dependent
- Complex symmetry handling
- Requires validation

## Application Areas
- Symmetry-constrained generation
- Crystal structure prediction
- Materials discovery
- Space group targeting

## Comparison with Other Codes
- **vs DiffCSP**: SyMat more symmetry-focused
- **vs CDVAE**: Different symmetry approach
- **Unique strength**: Explicit symmetry awareness, NeurIPS 2023

## Community and Support
- Academic development
- NeurIPS publication
- Research code

## Verification & Sources
**Primary sources**:
1. Paper: NeurIPS 2023
2. arXiv: https://arxiv.org/abs/2307.02707

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Documentation: AVAILABLE (paper)
- Development: ACADEMIC
- Applications: Symmetry-aware crystal generation
