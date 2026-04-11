# Skala

## Official Resources
- Source Repository: https://github.com/microsoft/skala
- License: Open source (MIT)

## Overview
**Skala** is a deep learning exchange-correlation functional developed by Microsoft that achieves chemical accuracy for atomization energies. It uses message-passing for non-local correlation and integrates with PySCF for DFT calculations.

**Scientific domain**: Deep learning XC functional with chemical accuracy  
**Target user community**: Researchers needing ML-enhanced DFT with chemical accuracy

## Theoretical Methods
- Deep learning XC functional
- Message-passing for non-local correlation
- Chemical accuracy (<1 kcal/mol)
- PySCF integration
- GPU and CPU implementations

## Capabilities (CRITICAL)
- Chemical accuracy atomization energies
- Non-local correlation via message passing
- PySCF integration
- GPU (Skala-GPU) and CPU (Skala-CPU)
- Broad thermochemistry benchmarks

**Sources**: GitHub repository, arXiv:2506.14665

## Key Strengths

### Accuracy:
- Chemical accuracy (<1 kcal/mol)
- Better than hybrid functionals
- Non-local correlation
- Broad benchmark performance

### Integration:
- PySCF self-consistent
- GPU acceleration
- CPU fallback
- Standard DFT workflow

## Inputs & Outputs
- **Input formats**: Molecular/periodic structures
- **Output data types**: XC energy, potential, total energy

## Interfaces & Ecosystem
- **PySCF**: DFT engine
- **JAX/PyTorch**: Backend
- **Python**: Core

## Performance Characteristics
- **Speed**: Similar to hybrid DFT
- **Accuracy**: Chemical accuracy
- **System size**: Molecular
- **Automation**: Full

## Computational Cost
- **SCF**: Similar to hybrid DFT
- **Training**: Days on GPU

## Limitations & Known Constraints
- **Molecular focus**: Primarily molecules
- **PySCF only**: Limited DFT code integration
- **New project**: Still maturing
- **GPU preferred**: For speed

## Comparison with Other Codes
- **vs DM21**: Skala has message-passing, DM21 is piecewise linear
- **vs DeePKS**: Skala is standalone, DeePKS is perturbative
- **vs NeuralXC**: Skala is more accurate, NeuralXC is simpler
- **Unique strength**: Deep learning XC functional achieving chemical accuracy with message-passing non-local correlation

## Application Areas

### Accurate DFT:
- Thermochemistry
- Kinetics
- Reaction energies
- Benchmark studies

### ML-DFT Research:
- XC functional development
- DFT accuracy improvement
- Method comparison

## Best Practices
- Use PySCF for SCF calculations
- Compare with standard functionals
- Validate on target chemistry
- Use GPU for speed

## Community and Support
- Open source (MIT)
- Microsoft maintained
- Azure AI integration

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/microsoft/skala

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Deep learning XC functional achieving chemical accuracy with message-passing non-local correlation
