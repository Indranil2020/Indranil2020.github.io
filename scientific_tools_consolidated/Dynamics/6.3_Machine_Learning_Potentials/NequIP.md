# NequIP

## Official Resources
- Homepage: https://github.com/mir-group/nequip
- Documentation: https://github.com/mir-group/nequip
- Source Repository: https://github.com/mir-group/nequip
- License: MIT

## Overview
NequIP (Neural Equivariant Interatomic Potentials) is an open-source code for building E(3)-equivariant interatomic potentials using graph neural networks. It achieves state-of-the-art accuracy with remarkable data efficiency, often requiring 100-1000x less training data than other ML potentials.

**Scientific domain**: Machine learning potentials, equivariant neural networks  
**Target user community**: Researchers developing accurate ML potentials with limited data

## Theoretical Methods
- E(3)-equivariant graph neural networks
- Equivariant convolutions
- Tensor product operations
- Message passing neural networks
- Invariant/equivariant features

## Capabilities (CRITICAL)
- E(3)-equivariant architecture
- Exceptional data efficiency
- High accuracy
- ASE calculator interface
- LAMMPS integration
- Multi-GPU training
- Flexible architecture

## Key Strengths

### Data Efficiency:
- 100-1000x less data needed
- Excellent for small datasets
- Strong extrapolation
- Few-shot learning capable

### Equivariance:
- Full E(3) equivariance
- Physical symmetries built-in
- Smooth energy surfaces
- Accurate forces

## Inputs & Outputs
- **Input formats**:
  - ASE Atoms
  - Extended XYZ
  - VASP OUTCAR
  
- **Output data types**:
  - Energies
  - Forces
  - Model checkpoints

## Interfaces & Ecosystem
- **ASE**: Calculator
- **LAMMPS**: pair_nequip
- **Allegro**: Scalable extension
- **PyTorch**: Backend

## Advanced Features
- **Equivariant convolutions**: E(3) symmetry
- **Data efficiency**: Few training points
- **Allegro extension**: Large-scale systems
- **Uncertainty**: Ensemble methods
- **Active learning**: Data selection

## Performance Characteristics
- Moderate inference speed
- Excellent accuracy/data tradeoff
- GPU acceleration
- Scales with model size

## Computational Cost
- Training: Hours (small datasets)
- Inference: Moderate speed
- Allegro faster for large systems
- Overall: Good for accuracy-critical applications

## Best Practices
- Start with small models
- Use data augmentation
- Validate carefully
- Consider Allegro for large systems

## Limitations & Known Constraints
- Slower than some alternatives
- GPU recommended
- Complex architecture
- Allegro better for large systems

## Application Areas
- Small molecule dynamics
- Materials with limited data
- High-accuracy requirements
- Method development
- Active learning

## Community and Support
- Active development (Harvard)
- GitHub issues
- Growing community
- Good documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mir-group/nequip
2. S. Batzner et al., Nat. Commun. 13, 2453 (2022)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Published in Nature Communications
