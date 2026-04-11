# MatterSim

## Official Resources
- Source Repository: https://github.com/microsoft/mattersim
- License: Open source (MIT)

## Overview
**MatterSim** is Microsoft's deep learning interatomic potential trained across elements, temperatures, and pressures. It covers 94 elements with multi-domain training data, achieving state-of-art accuracy on diverse materials benchmarks.

**Scientific domain**: Universal deep learning potential across T, P, and compositions  
**Target user community**: Researchers needing universal potential for diverse conditions

## Theoretical Methods
- Graph neural network architecture
- Multi-domain training (T, P, composition)
- Foundation model approach
- 94 elements coverage
- Fine-tuning capability

## Capabilities (CRITICAL)
- Universal potential (94 elements)
- Temperature and pressure coverage
- Fine-tuning for specific systems
- ASE calculator
- LAMMPS integration
- Matbench Discovery benchmarked

**Sources**: GitHub repository, arXiv:2405.04967

## Key Strengths

### Multi-Domain:
- Temperature coverage
- Pressure coverage
- Diverse compositions
- Phase transitions

### Universal:
- 94 elements
- Pretrained foundation model
- No retraining needed
- Fine-tuning available

### Microsoft:
- Well-resourced development
- Regular updates
- Azure integration
- Professional support

## Inputs & Outputs
- **Input formats**: Structures (ASE/pymatgen)
- **Output data types**: Energies, forces, stresses

## Interfaces & Ecosystem
- **ASE**: Calculator
- **LAMMPS**: MD engine
- **PyTorch**: Backend

## Performance Characteristics
- **Speed**: Fast (GPU)
- **Accuracy**: Near-DFT across conditions
- **System size**: 1-10000+ atoms
- **Automation**: Full

## Computational Cost
- **MD**: ~1000x faster than DFT
- **Fine-tuning**: Minutes to hours

## Limitations & Known Constraints
- **PBE-level**: Trained on PBE data
- **GPU required**: For production
- **New project**: Still maturing
- **Microsoft license**: Some restrictions

## Comparison with Other Codes
- **vs MACE-MP-0**: MatterSim has T/P coverage, MACE is room T
- **vs CHGNet**: MatterSim is multi-domain, CHGNet is charge-aware
- **vs ORB**: MatterSim is Microsoft, ORB is independent
- **Unique strength**: Multi-domain universal potential covering temperature, pressure, and 94 elements

## Application Areas

### Multi-Condition MD:
- High-temperature MD
- High-pressure simulations
- Phase diagrams
- Melting and crystallization

### Universal:
- General-purpose MD
- Structure relaxation
- Energy screening

## Best Practices
- Use pretrained model first
- Fine-tune for specific conditions
- Validate at target T/P
- Compare with DFT MD

## Community and Support
- Open source (MIT)
- Microsoft maintained
- GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/microsoft/mattersim

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Multi-domain universal potential covering temperature, pressure, and 94 elements
