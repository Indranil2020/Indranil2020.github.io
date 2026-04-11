# FAIR-Chem (Open Catalyst Project / OCP)

## Official Resources
- Source Repository: https://github.com/FAIR-Chem/fairchem
- Documentation: https://fair-chem.github.io/
- PyPI: https://pypi.org/project/fairchem-core/
- License: Open source (MIT)

## Overview
**FAIR-Chem** (formerly Open Catalyst Project) is Meta's AI research initiative for catalysis and materials science. It provides pretrained equivariant models (EquiformerV2, eSCN, GemNet-OC) trained on OC20/OC22 datasets, with a comprehensive Python package for training and inference.

**Scientific domain**: Catalysis ML, equivariant models, large-scale datasets  
**Target user community**: Researchers in catalysis, surface science, and ML for chemistry

## Theoretical Methods
- EquiformerV2 (equivariant transformer)
- eSCN (equivariant spherical channel network)
- GemNet-OC (geometric message passing)
- DimeNet++ (directional message passing)
- SchNet (continuous-filter convolution)
- PaiNN (polarizable atom interaction)

## Capabilities (CRITICAL)
- Pretrained models for catalysis (OC20, OC22)
- Training framework for custom models
- Multiple architecture support
- ASE calculator interface
- LAMMPS integration
- 130K+ monthly PyPI downloads

**Sources**: GitHub repository, ACS Catal. 11, 6059 (2021)

## Key Strengths

### Catalysis Focus:
- OC20: 1.3M+ DFT relaxations
- OC22: Oxide surfaces
- Adsorption energy prediction
- Reaction pathway modeling

### Multiple Architectures:
- EquiformerV2 (best accuracy)
- eSCN (efficient)
- GemNet-OC (balanced)
- SchNet, PaiNN (baselines)

### Production Ready:
- PyPI installable
- 130K+ monthly downloads
- Active development (69 contributors)
- Comprehensive documentation

## Inputs & Outputs
- **Input formats**: Atomic structures, datasets
- **Output data types**: Energies, forces, adsorption energies

## Interfaces & Ecosystem
- **ASE**: Calculator
- **LAMMPS**: MD engine
- **PyTorch**: Backend
- **Python**: Core

## Performance Characteristics
- **Speed**: Fast (GPU)
- **Accuracy**: State-of-art on OC20/OC22
- **System size**: Surface slabs
- **Automation**: Full

## Computational Cost
- **Inference**: Milliseconds
- **Training**: Days on multi-GPU
- **Pretrained**: Available

## Limitations & Known Constraints
- **Catalysis focus**: Optimized for surfaces
- **GPU required**: For training
- **Large package**: Many dependencies
- **Dataset specific**: OC20/OC22 benchmarks

## Comparison with Other Codes
- **vs MACE**: FAIR-Chem is catalysis, MACE is universal
- **vs SchNetPack**: FAIR-Chem has pretrained models, SchNetPack is framework
- **vs TorchANI**: FAIR-Chem is periodic, TorchANI is molecular
- **Unique strength**: Pretrained catalysis models on OC20/OC22 datasets with multiple architectures

## Application Areas

### Catalysis:
- Adsorption energy prediction
- Catalyst screening
- Surface reactions
- OC20/OC22 benchmarking

### Model Development:
- Custom architecture training
- Transfer learning
- Active learning for catalysis

## Best Practices

### Usage:
- Use pretrained EquiformerV2 for best accuracy
- Fine-tune on target systems
- Use OC20/OC22 for benchmarking
- Validate with DFT

## Community and Support
- Open source (MIT)
- PyPI installable (fairchem-core)
- Meta AI maintained
- 69 contributors
- Comprehensive documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/FAIR-Chem/fairchem

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE (130K+ monthly downloads)
- Specialized strength: Pretrained catalysis models on OC20/OC22 with multiple equivariant architectures
