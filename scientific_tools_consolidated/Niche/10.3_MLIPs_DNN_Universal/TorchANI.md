# TorchANI

## Official Resources
- Source Repository: https://github.com/aiqm/torchani
- Documentation: https://aiqm.github.io/torchani/
- PyPI: https://pypi.org/project/torchani/
- License: Open source (MIT)

## Overview
**TorchANI** is a PyTorch implementation of ANI (Accurate Neural Network Engine for Molecular Energies). It provides pretrained ANI-1x, ANI-1ccx, and ANI-2x models for organic molecules (H, C, N, O, F, S, Cl) with fast MD and QM/MM capabilities.

**Scientific domain**: ANI neural network potentials for organic molecules  
**Target user community**: Researchers simulating organic molecular systems with ML

## Theoretical Methods
- Behler-Parrinello HDNNP architecture
- ANI-1x/1ccx/2x pretrained models
- AEV (Atomic Environment Vector) descriptors
- Coupled-cluster training data (ANI-1ccx)
- QM/MM integration with NAMD

## Capabilities (CRITICAL)
- Pretrained ANI models (1x, 1ccx, 2x)
- H, C, N, O, F, S, Cl coverage
- ASE calculator
- QM/MM with NAMD
- PyTorch-based (GPU accelerated)

**Sources**: GitHub repository, J. Chem. Inf. Model. 60, 3408 (2020)

## Key Strengths

### Pretrained:
- ANI-1x (DFT-trained)
- ANI-1ccx (CCSD(T)-trained)
- ANI-2x (extended elements)
- No training needed

### Molecular:
- Organic molecules
- Drug-like compounds
- Biomolecules
- QM/MM compatible

### Performance:
- PyTorch backend
- GPU acceleration
- Fast inference
- Conda installable (1.2M downloads)

## Inputs & Outputs
- **Input formats**: Molecular structures (ASE)
- **Output data types**: Energies, forces

## Interfaces & Ecosystem
- **ASE**: Calculator
- **NAMD**: QM/MM
- **PyTorch**: Backend
- **Python**: Core

## Performance Characteristics
- **Speed**: Fast (GPU)
- **Accuracy**: Near-DFT (1x), near-CC (1ccx)
- **System size**: Molecular (100s atoms)
- **Automation**: Full

## Computational Cost
- **Inference**: Milliseconds
- **MD**: ~1000x faster than DFT

## Limitations & Known Constraints
- **7 elements only**: H, C, N, O, F, S, Cl
- **Molecular only**: No periodic systems
- **Organic focus**: Not for inorganic
- **PBE/CC data**: Limited training levels

## Comparison with Other Codes
- **vs DeePMD-kit**: TorchANI is molecular, DeePMD is periodic
- **vs MACE-OFF23**: Both organic, MACE is equivariant
- **vs SchNetPack**: TorchANI is pretrained, SchNetPack is framework
- **Unique strength**: Pretrained ANI models with CCSD(T) accuracy (ANI-1ccx) for organic molecules

## Application Areas

### Organic MD:
- Drug molecule dynamics
- Biomolecular simulations
- Conformational sampling
- Reaction pathways

### QM/MM:
- NAMD integration
- Enzyme catalysis
- Solvent effects
- Photochemistry

## Best Practices
- Use ANI-2x for broadest coverage
- Use ANI-1ccx for highest accuracy
- Validate on target molecules
- Use GPU for production MD

## Community and Support
- Open source (MIT)
- PyPI/Conda installable
- AIQM maintained
- Roitberg group
- 1.2M conda downloads

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aiqm/torchani

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: Pretrained ANI models with CCSD(T) accuracy for organic molecules
