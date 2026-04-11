# DeePKS-kit

## Official Resources
- Source Repository: https://github.com/deepmodeling/deepks-kit
- License: Open source (LGPL-3.0)

## Overview
**DeePKS-kit** is a package for developing machine-learned XC functionals for quantum chemistry. It implements both perturbative (DeePHF) and self-consistent (DeePKS) schemes, integrating with PySCF for accurate energy functionals.

**Scientific domain**: ML XC functional development (perturbative and self-consistent)  
**Target user community**: Researchers developing and applying ML XC functionals

## Theoretical Methods
- DeePKS (self-consistent scheme)
- DeePHF (perturbative scheme)
- Neural network XC functional
- PySCF integration
- ABACUS integration

## Capabilities (CRITICAL)
- Self-consistent DeePKS
- Perturbative DeePHF
- PySCF and ABACUS integration
- Neural network XC training
- Transferable functionals

**Sources**: GitHub repository

## Key Strengths

### Dual Scheme:
- Self-consistent (DeePKS)
- Perturbative (DeePHF)
- Both approaches available
- Accuracy vs cost tradeoff

### Integration:
- PySCF (molecular)
- ABACUS (periodic)
- DP-GEN workflow
- DeepModeling ecosystem

## Inputs & Outputs
- **Input formats**: Training data (QM calculations)
- **Output data types**: Trained XC models, corrected energies

## Interfaces & Ecosystem
- **PySCF**: Molecular DFT
- **ABACUS**: Periodic DFT
- **DP-GEN**: Active learning
- **Python**: Core

## Performance Characteristics
- **Speed**: Similar to hybrid DFT (SC)
- **Accuracy**: Near coupled-cluster
- **System size**: Molecular and periodic
- **Automation**: Full

## Computational Cost
- **SC calculation**: Similar to hybrid
- **Training**: Hours on GPU

## Limitations & Known Constraints
- **Training required**: Need reference data
- **PySCF/ABACUS only**: Limited code support
- **Complex setup**: Multi-step workflow
- **LGPL license**: Copyleft

## Comparison with Other Codes
- **vs Skala**: DeePKS is perturbative+SC, Skala is message-passing
- **vs DM21**: DeePKS is trainable, DM21 is fixed
- **vs NeuralXC**: DeePKS is more sophisticated
- **Unique strength**: Dual perturbative/self-consistent ML XC with PySCF and ABACUS integration

## Application Areas

### Accurate DFT:
- Beyond-DFT accuracy at DFT cost
- Molecular property prediction
- Periodic system calculations
- Multi-level calculations

### XC Development:
- Custom XC functional training
- Active learning for XC
- Transfer learning

## Best Practices
- Start with perturbative scheme
- Validate against high-level QM
- Use DP-GEN for training data
- Fine-tune for target chemistry

## Community and Support
- Open source (LGPL-3.0)
- DeepModeling community
- GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/deepmodeling/deepks-kit

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Dual perturbative/self-consistent ML XC with PySCF and ABACUS integration
