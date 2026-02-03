# OMatG

## Overview
OMatG (Open Materials Generation) is a state-of-the-art generative model for crystal structure prediction and de novo generation of inorganic crystals using stochastic interpolants.

## Theoretical Basis
- Stochastic interpolants
- Flow-based generation
- Crystal structure prediction
- De novo generation
- Equivariant architecture

## Key Capabilities
- Crystal structure prediction
- De novo crystal generation
- State-of-the-art performance
- Flexible generation modes
- ICML 2025 publication

**Sources**: ICML 2025, OpenReview

## Key Strengths

### Methodology:
- Stochastic interpolants
- Flexible framework
- Multiple generation modes

### Performance:
- State-of-the-art results
- Outperforms flow/diffusion baselines
- Well-benchmarked

### Flexibility:
- CSP mode
- De novo mode
- Adaptable architecture

## Inputs & Outputs
- **Input formats**: Composition (CSP), nothing (de novo)
- **Output data types**: Generated crystal structures

## Interfaces & Ecosystem
- **Framework**: PyTorch Lightning
- **Hugging Face**: Model checkpoints
- **Datasets**: Benchmark datasets

## Workflow and Usage
1. Load pretrained model
2. Select generation mode
3. Run generation
4. Evaluate structures
5. Validate with DFT

## Performance Characteristics
- GPU-accelerated
- Efficient generation
- High quality outputs

## Computational Cost
- Generation: fast
- Training: GPU hours
- Validation: DFT

## Best Practices
- Use appropriate generation mode
- Generate multiple samples
- Validate predictions
- Check structural validity

## Limitations & Known Constraints
- Training data dependent
- Recent development
- Requires validation

## Application Areas
- Crystal structure prediction
- De novo material generation
- Materials discovery
- Generative materials science

## Comparison with Other Codes
- **vs FlowMM**: OMatG stochastic interpolants
- **vs DiffCSP**: Different generative approach
- **Unique strength**: Stochastic interpolants, ICML 2025

## Community and Support
- Open-source (GitHub)
- Hugging Face models
- Academic development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/FERMat-ML/OMatG
2. Hugging Face: https://huggingface.co/OMatG
3. Paper: ICML 2025

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: ACTIVE
- Applications: Crystal structure generation
