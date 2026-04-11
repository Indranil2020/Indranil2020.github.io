# SevenNet

## Official Resources
- Source Repository: https://github.com/MDIL-SNU/SevenNet
- Documentation: https://sevennet.readthedocs.io/
- PyPI: https://pypi.org/project/sevennet/
- License: Open source (GPL-3.0)

## Overview
**SevenNet** (Scalable Equivariant NEural network potential) is a graph neural network interatomic potential with equivariant architecture. It provides pretrained universal models (SevenNet-0) trained on the MPF/MPtrj dataset, supporting fine-tuning and LAMMPS integration.

**Scientific domain**: Equivariant GNN universal interatomic potential  
**Target user community**: Researchers needing scalable equivariant MLIP for MD

## Theoretical Methods
- E(3)-equivariant message passing
- Scalable architecture
- Pretrained universal model (SevenNet-0)
- Fine-tuning capability
- Parallel MD support

## Capabilities (CRITICAL)
- Pretrained universal potential (89+ elements)
- Fine-tuning interface
- ASE calculator
- LAMMPS integration (pair_style)
- Parallel GNN-MD simulation
- Matbench Discovery benchmarked

**Sources**: GitHub repository, arXiv:2409.04649

## Key Strengths

### Equivariant:
- E(3)-equivariance
- Rotational invariance
- Better data efficiency
- Physically consistent

### Scalable:
- Parallel MD support
- GPU acceleration
- Large system support
- LAMMPS integration

### Pretrained:
- SevenNet-0 universal model
- MPF/MPtrj training data
- No training needed
- Fine-tuning available

## Inputs & Outputs
- **Input formats**: Structures (ASE/pymatgen)
- **Output data types**: Energies, forces, stresses

## Interfaces & Ecosystem
- **ASE**: Calculator
- **LAMMPS**: MD engine
- **PyTorch**: Backend
- **Python**: Core

## Performance Characteristics
- **Speed**: Fast (GPU)
- **Accuracy**: Near-DFT
- **System size**: 1-100000+ atoms (parallel)
- **Automation**: Full

## Computational Cost
- **MD**: ~1000x faster than DFT
- **Fine-tuning**: Minutes to hours

## Limitations & Known Constraints
- **GPL license**: Copyleft
- **PBE-level**: Trained on PBE data
- **GPU required**: For production
- **New project**: Still maturing

## Comparison with Other Codes
- **vs MACE**: SevenNet is scalable, MACE is more accurate
- **vs NequIP**: SevenNet is pretrained, NequIP needs training
- **vs CHGNet**: SevenNet is equivariant, CHGNet has charge
- **Unique strength**: Scalable equivariant universal potential with parallel MD support

## Application Areas

### Large-Scale MD:
- Parallel molecular dynamics
- High-throughput screening
- Phase transitions
- Mechanical properties

### Fine-Tuning:
- Element-specific accuracy
- New chemistry coverage
- Active learning

## Best Practices

### Usage:
- Start with SevenNet-0
- Fine-tune for specific systems
- Use LAMMPS for production MD
- Validate against DFT

## Community and Support
- Open source (GPL-3.0)
- PyPI installable
- MDIL-SNU maintained
- ReadTheDocs documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/MDIL-SNU/SevenNet

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: Scalable equivariant universal potential with parallel MD and fine-tuning
