# CrystalFlow

## Overview
CrystalFlow is a flow-based generative model for crystalline materials using Continuous Normalizing Flows (CNFs) within the Conditional Flow Matching (CFM) framework. Published in Nature Communications 2025.

## Theoretical Basis
- Continuous Normalizing Flows (CNFs)
- Conditional Flow Matching (CFM)
- Periodic structure handling
- Lattice and coordinate generation
- Property-guided generation

## Key Capabilities
- Crystal structure generation
- Flow-based sampling
- Property optimization
- Stable structure generation
- Nature Communications publication

**Sources**: Nature Communications (2025), arXiv:2412.11693

## Key Strengths

### Methodology:
- Flow-based (not diffusion)
- CNF framework
- Efficient sampling

### Performance:
- Competitive results
- Multiple benchmarks
- Good stability

### Publication:
- Nature Communications 2025
- Peer-reviewed
- Well-documented

## Inputs & Outputs
- **Input formats**: Training structures, composition (generation)
- **Output data types**: Generated crystal structures

## Interfaces & Ecosystem
- **Framework**: PyTorch
- **Datasets**: Standard benchmarks
- **Evaluation**: DFT validation

## Workflow and Usage
1. Train on crystal dataset
2. Sample from flow model
3. Generate structures
4. Evaluate stability
5. Validate with DFT

## Performance Characteristics
- GPU-accelerated
- Efficient flow sampling
- Good generation quality

## Computational Cost
- Training: GPU hours
- Sampling: fast
- Validation: DFT

## Best Practices
- Use appropriate training data
- Generate multiple samples
- Validate predictions
- Check structural validity

## Limitations & Known Constraints
- Training data dependent
- Recent development
- Requires validation

## Application Areas
- Crystal structure generation
- Materials discovery
- Property-guided design
- Generative materials science

## Comparison with Other Codes
- **vs FlowMM**: Similar approach, different implementation
- **vs DiffCSP**: Flow vs diffusion
- **vs CDVAE**: Different generative framework
- **Unique strength**: CNF/CFM framework, Nature Comms 2025

## Community and Support
- Open-source (GitHub)
- Academic development
- Recent publication

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ixsluo/CrystalFlow
2. Publication: Nature Communications (2025)
3. arXiv: https://arxiv.org/abs/2412.11693

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: RECENT (2025)
- Applications: Crystal structure generation
