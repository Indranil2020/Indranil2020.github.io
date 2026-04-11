# ORB

## Official Resources
- Source Repository: https://github.com/orbital-materials/orb-models
- License: Open source (Apache-2.0)

## Overview
**ORB** (Open Reusable Bindings) is an open-source universal interatomic potential covering 117 elements with 25M parameters. Trained on multiple datasets, it achieves competitive accuracy on Matbench Discovery with efficient inference.

**Scientific domain**: Universal potential with broadest element coverage (117)  
**Target user community**: Researchers needing widest element coverage in a universal potential

## Theoretical Methods
- Graph neural network architecture
- 117 elements coverage
- 25M parameters
- Multi-dataset training
- Efficient inference

## Capabilities (CRITICAL)
- 117 elements (broadest coverage)
- 25M parameter model
- Matbench Discovery benchmarked
- ASE calculator
- Efficient GPU inference

**Sources**: GitHub repository, arXiv:2502.20851

## Key Strengths

### Broadest Coverage:
- 117 elements
- Most of periodic table
- Rare earth elements
- Actinides

### Performance:
- 25M parameters
- Competitive accuracy
- Efficient inference
- Matbench Discovery tested

## Inputs & Outputs
- **Input formats**: Structures (ASE)
- **Output data types**: Energies, forces, stresses

## Interfaces & Ecosystem
- **ASE**: Calculator
- **PyTorch**: Backend

## Performance Characteristics
- **Speed**: Fast (GPU)
- **Accuracy**: Competitive
- **System size**: 1-10000+ atoms

## Computational Cost
- **MD**: ~1000x faster than DFT
- **Inference**: Milliseconds

## Limitations & Known Constraints
- **New project**: Still maturing
- **GPU required**: For production
- **Limited MD integration**: ASE only
- **Documentation**: Growing

## Comparison with Other Codes
- **vs MACE-MP-0**: ORB covers 117 elements, MACE covers 89
- **vs CHGNet**: ORB has broader coverage, CHGNet has charge
- **Unique strength**: Broadest element coverage (117) with 25M parameter universal potential

## Application Areas

### Wide-Coverage MD:
- Rare earth materials
- Actinide simulations
- Multi-element alloys
- Unexplored chemistries

### Screening:
- High-throughput energy evaluation
- Structure relaxation
- Stability prediction

## Best Practices
- Use for systems with rare elements
- Validate against DFT
- Compare with other UIPs

## Community and Support
- Open source (Apache-2.0)
- Orbital Materials maintained
- GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/orbital-materials/orb-models

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Broadest element coverage (117) with 25M parameter universal potential
