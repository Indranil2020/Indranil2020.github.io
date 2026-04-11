# PaiNN

## Official Resources
- Source Repository: https://github.com/atomistic-machine-learning/schnetpack (PaiNN module)
- Paper: ICML 2021
- License: Open source (MIT)

## Overview
**PaiNN** (Polarizable Atom Interaction Neural Network) is an equivariant message passing architecture that uses scalar and vector features. It achieves high accuracy for molecular property prediction with efficient equivariant message passing.

**Scientific domain**: Equivariant message passing with polarizable atom interactions  
**Target user community**: Researchers needing efficient equivariant model for molecular properties

## Theoretical Methods
- Equivariant message passing
- Scalar and vector features
- Polarizable atom interactions
- Efficient equivariance (no spherical harmonics)
- SchNetPack integration

## Capabilities (CRITICAL)
- Equivariant predictions
- Scalar + vector features
- Molecular property prediction
- SchNetPack integration
- Efficient architecture

**Sources**: SchNetPack repository, ICML 2021

## Key Strengths

### Efficient Equivariance:
- No spherical harmonics
- Simple vector features
- Fast training
- Good data efficiency

### Molecular:
- QM9 benchmarks
- Molecular dynamics
- Property prediction
- SchNetPack ecosystem

## Inputs & Outputs
- **Input formats**: Molecular structures
- **Output data types**: Energies, forces, dipole moments, etc.

## Interfaces & Ecosystem
- **SchNetPack**: Framework
- **ASE**: Calculator
- **PyTorch**: Backend

## Performance Characteristics
- **Speed**: Fast (efficient equivariance)
- **Accuracy**: State-of-art on QM9
- **System size**: Molecular
- **Automation**: Full

## Computational Cost
- **Training**: Hours on GPU
- **Inference**: Milliseconds

## Limitations & Known Constraints
- **Molecular focus**: Primarily non-periodic
- **SchNetPack only**: No standalone
- **Not universal**: Needs training

## Comparison with Other Codes
- **vs NequIP**: PaiNN is simpler, NequIP is higher-order
- **vs SchNet**: PaiNN is equivariant, SchNet is invariant
- **Unique strength**: Efficient equivariant message passing without spherical harmonics

## Application Areas

### Molecular Properties:
- QM9 benchmark
- Dipole moments
- Polarizabilities
- HOMO-LUMO gaps

### MD:
- Molecular dynamics
- Conformational sampling
- Free energy calculations

## Best Practices
- Use SchNetPack for training
- Start with QM9 for benchmarking
- Fine-tune for target chemistry

## Community and Support
- Open source (MIT)
- SchNetPack ecosystem
- ICML published

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/atomistic-machine-learning/schnetpack

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (via SchNetPack)
- Specialized strength: Efficient equivariant message passing without spherical harmonics
