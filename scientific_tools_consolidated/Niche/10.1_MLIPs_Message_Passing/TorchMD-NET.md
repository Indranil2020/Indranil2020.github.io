# TorchMD-NET

## Official Resources
- Source Repository: https://github.com/torchmd/torchmd-net
- Documentation: https://torchmd-net.readthedocs.io/
- Conda: https://anaconda.org/conda-forge/torchmd-net
- License: Open source (MIT)

## Overview
**TorchMD-NET** is a PyTorch-based neural network potential for molecular dynamics. It implements multiple architectures including TorchMD-Net, PaiNN, and Transformer-based models, with pretrained models and efficient MD integration.

**Scientific domain**: Neural network potentials for molecular dynamics  
**Target user community**: Researchers running ML-driven molecular dynamics

## Theoretical Methods
- TorchMD-Net architecture
- PaiNN (polarizable atom interaction)
- Transformer-based potentials
- Tensor field networks

## Capabilities (CRITICAL)
- Multiple architecture support
- Pretrained models
- TorchMD integration for MD
- ASE calculator
- Conda installable (790K downloads)

**Sources**: GitHub repository, arXiv:2202.02541

## Key Strengths

### Multi-Architecture:
- PaiNN, TorchMD-Net, Transformers
- Easy architecture switching
- Benchmarking framework

### MD Integration:
- TorchMD for simulations
- ASE calculator
- Efficient inference

### Popular:
- 790K+ conda downloads
- Active development
- Well documented

## Inputs & Outputs
- **Input formats**: Molecular structures
- **Output data types**: Energies, forces

## Interfaces & Ecosystem
- **TorchMD**: MD engine
- **ASE**: Calculator
- **PyTorch**: Backend

## Performance Characteristics
- **Speed**: Fast (GPU)
- **Accuracy**: Architecture-dependent
- **System size**: Molecular
- **Automation**: Full

## Computational Cost
- **Inference**: Milliseconds
- **Training**: Hours on GPU

## Limitations & Known Constraints
- **Molecular focus**: Primarily non-periodic
- **No universal model**: Need training
- **TorchMD dependency**: For MD

## Comparison with Other Codes
- **vs SchNetPack**: TorchMD-NET is MD-focused, SchNetPack is broader
- **vs MACE**: TorchMD-NET is molecular, MACE is universal
- **Unique strength**: Multi-architecture NN potential with integrated MD (TorchMD)

## Application Areas

### Molecular MD:
- ML-driven molecular dynamics
- Property prediction
- Conformational sampling

### Architecture Research:
- Benchmarking architectures
- PaiNN vs Transformer
- Model comparison

## Best Practices
- Start with PaiNN for balanced performance
- Use pretrained models when available
- Validate with QM benchmarks

## Community and Support
- Open source (MIT)
- Conda installable
- TorchMD team maintained
- ReadTheDocs documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/torchmd/torchmd-net

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Conda: AVAILABLE (790K downloads)
- Specialized strength: Multi-architecture NN potential with integrated MD (TorchMD)
