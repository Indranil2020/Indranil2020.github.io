# CHGNet

## Official Resources
- Source Repository: https://github.com/CederGroupHub/chgnet
- Documentation: https://chgnet.readthedocs.io/
- PyPI: https://pypi.org/project/chgnet/
- License: Open source (MIT)

## Overview
**CHGNet** (Crystal Hamiltonian Graph neural Network) is a pretrained universal neural network potential that incorporates charge information. It achieves near-DFT accuracy for materials simulations with explicit electron interactions and charge distribution modeling.

**Scientific domain**: Universal MLIP with charge equilibration for materials  
**Target user community**: Researchers needing fast, accurate MD with charge-aware potentials

## Theoretical Methods
- Graph neural network with charge-informed features
- Universal potential covering 94 elements
- Magnetic moment prediction
- Charge distribution modeling
- Materials Project training data

## Capabilities (CRITICAL)
- Universal potential for 94 elements
- Charge-aware predictions
- Magnetic moment prediction
- ASE calculator interface
- LAMMPS integration
- Matbench Discovery benchmarked

**Sources**: GitHub repository, arXiv:2503.09814

## Key Strengths

### Charge-Aware:
- Explicit charge modeling
- Electron interactions
- Better redox chemistry
- Ionic bonding accuracy

### Universal:
- 94 elements covered
- Pretrained on MP data
- No retraining needed for most systems
- Fine-tuning supported

### Integration:
- ASE calculator
- LAMMPS pair_style
- pymatgen compatibility
- Matbench Discovery

## Inputs & Outputs
- **Input formats**: Structures (pymatgen/ASE)
- **Output data types**: Energies, forces, stresses, magnetic moments, charges

## Interfaces & Ecosystem
- **ASE**: Calculator interface
- **LAMMPS**: MD engine
- **pymatgen**: Structure handling
- **Python/PyTorch**: Core

## Performance Characteristics
- **Speed**: ~1ms/atom (GPU)
- **Accuracy**: Near-DFT (MAE ~25 meV/atom)
- **System size**: 1-10000+ atoms
- **Automation**: Full

## Computational Cost
- **MD**: ~1000x faster than DFT
- **Training**: Hours on GPU (pretrained available)

## Limitations & Known Constraints
- **PBE-level accuracy**: Trained on PBE data
- **Charge model limitations**: Approximate
- **GPU preferred**: CPU is slower
- **Not for molecules**: Optimized for periodic

## Comparison with Other Codes
- **vs MACE-MP-0**: CHGNet has charge, MACE has higher accuracy
- **vs M3GNet**: CHGNet adds charge, M3GNet is simpler
- **vs GAP**: CHGNet is universal, GAP is element-specific
- **Unique strength**: Charge-aware universal potential with magnetic moment prediction

## Application Areas

### Materials MD:
- High-throughput MD simulations
- Phase stability prediction
- Ionic conductivity
- Redox chemistry

### Structure Prediction:
- Relaxation from unrelaxed structures
- Energy above hull prediction
- Metastable phase exploration

## Best Practices

### Usage:
- Use pretrained model first
- Fine-tune for specific chemistries
- Validate against DFT for critical systems
- Use GPU for production MD

## Community and Support
- Open source (MIT)
- PyPI installable
- Ceder Group maintained
- ReadTheDocs documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/CederGroupHub/chgnet

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: Charge-aware universal potential with magnetic moment prediction for 94 elements
