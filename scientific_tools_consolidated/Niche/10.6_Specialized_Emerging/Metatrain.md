# Metatrain

## Official Resources
- Source Repository: https://github.com/metatensor/metatrain
- Documentation: https://metatrain.readthedocs.io/
- PyPI: https://pypi.org/project/metatrain/
- License: Open source (BSD-3)

## Overview
**Metatrain** is a modular training framework for ML interatomic potentials. It provides a unified interface for training multiple architectures (MACE, PET, etc.) with the metatensor data format, supporting fine-tuning of foundation models like PET-MAD.

**Scientific domain**: Modular MLIP training with metatensor format  
**Target user community**: Researchers training and fine-tuning MLIPs across architectures

## Theoretical Methods
- Modular training framework
- Multiple architecture support (MACE, PET, etc.)
- Metatensor data format
- Foundation model fine-tuning
- Hyperparameter optimization

## Capabilities (CRITICAL)
- Train MACE, PET, and other architectures
- Fine-tune PET-MAD foundation model
- Metatensor data format
- Hyperparameter optimization
- 10K monthly PyPI downloads

**Sources**: GitHub repository

## Key Strengths

### Modular:
- Multiple architectures
- Same interface for all
- Easy architecture comparison
- Plugin system

### Fine-Tuning:
- PET-MAD foundation model
- Custom dataset adaptation
- Transfer learning
- Efficient training

### Metatensor:
- Standardized data format
- Efficient I/O
- Cross-architecture compatible
- Community standard

## Inputs & Outputs
- **Input formats**: Metatensor format, extended XYZ
- **Output data types**: Trained models, metatensor models

## Interfaces & Ecosystem
- **metatensor**: Data format
- **metatomic**: Model format
- **ASE**: Structure handling
- **Python**: Core

## Performance Characteristics
- **Speed**: Architecture-dependent
- **Accuracy**: Architecture-dependent
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Training**: Hours on GPU
- **Fine-tuning**: Minutes to hours

## Limitations & Known Constraints
- **New framework**: Still maturing
- **Limited architectures**: Growing list
- **Metatensor format**: Learning curve
- **Documentation**: Growing

## Comparison with Other Codes
- **vs MACE standalone**: Metatrain is multi-architecture
- **vs DP-GEN**: Metatrain is training, DP-GEN is active learning
- **vs SchNetPack**: Metatrain is modular, SchNetPack is monolithic
- **Unique strength**: Modular training framework for multiple MLIP architectures with metatensor format

## Application Areas

### MLIP Training:
- Multi-architecture comparison
- Foundation model fine-tuning
- Custom potential development
- Benchmark studies

### Research:
- Architecture development
- Transfer learning studies
- Data format standardization

## Best Practices
- Start with pretrained foundation model
- Fine-tune on target data
- Compare multiple architectures
- Use metatensor for data management

## Community and Support
- Open source (BSD-3)
- PyPI installable
- Metatensor community
- 34 contributors
- ReadTheDocs documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/metatensor/metatrain

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE (10K monthly)
- Specialized strength: Modular training framework for multiple MLIP architectures with metatensor format
