# EquiformerV2

## Official Resources
- Source Repository: https://github.com/atomicarchitects/equiformer_v2
- Paper: ICLR 2024
- License: Open source (MIT)

## Overview
**EquiformerV2** is an improved equivariant transformer for atomistic systems that scales to higher-degree representations. It achieves state-of-the-art performance on OC20 benchmarks and is the backbone model for FAIR-Chem's catalysis predictions.

**Scientific domain**: Equivariant transformer for interatomic potentials  
**Target user community**: Researchers needing highest-accuracy equivariant MLIPs

## Theoretical Methods
- Equivariant transformer architecture
- Higher-degree irreducible representations
- Attention mechanism with SO(3) equivariance
- Scalable to 153M parameters

## Capabilities (CRITICAL)
- State-of-art OC20 performance
- Higher-degree representations
- Pretrained models available
- FAIR-Chem integration
- Scalable architecture

**Sources**: GitHub repository, ICLR 2024

## Key Strengths

### Accuracy:
- Best on OC20 benchmarks
- Higher-degree equivariance
- Large model scaling (153M)

### Scalability:
- Multi-GPU training
- Large batch support
- Efficient attention

## Inputs & Outputs
- **Input formats**: Atomic structures
- **Output data types**: Energies, forces

## Interfaces & Ecosystem
- **FAIR-Chem**: Integration
- **PyTorch**: Backend
- **Python**: Core

## Performance Characteristics
- **Speed**: Moderate (large model)
- **Accuracy**: State-of-art
- **System size**: Surface slabs
- **Automation**: Full

## Computational Cost
- **Training**: Days on multi-GPU
- **Inference**: Seconds per structure

## Limitations & Known Constraints
- **Large model**: 153M parameters
- **GPU required**: For training and inference
- **Catalysis optimized**: OC20/OC22 focus

## Comparison with Other Codes
- **vs Equiformer**: V2 has higher-degree, better accuracy
- **vs MACE**: EquiformerV2 is catalysis, MACE is universal
- **Unique strength**: Highest accuracy equivariant transformer scaling to 153M parameters

## Application Areas

### Catalysis:
- OC20/OC22 benchmarking
- Adsorption energy prediction
- Surface reaction modeling

### MLIP Development:
- Architecture benchmarking
- Scaling studies
- Transfer learning

## Best Practices
- Use pretrained models
- Fine-tune for specific systems
- Multi-GPU for training

## Community and Support
- Open source (MIT)
- Atomic Architects maintained
- FAIR-Chem ecosystem

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/atomicarchitects/equiformer_v2

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Highest accuracy equivariant transformer scaling to 153M parameters
