# I-ReaxFF

## Official Resources
- Source Repository: https://github.com/fenggo/I-ReaxFF
- License: Open source

## Overview
**I-ReaxFF** (Intelligent-Reactive Force Field) is a differentiable ReaxFF framework based on TensorFlow. It enables gradient-based optimization and neural network augmentation of ReaxFF, combining reactive force fields with message passing neural networks.

**Scientific domain**: Differentiable ReaxFF with neural network augmentation  
**Target user community**: Researchers combining ReaxFF with neural networks for reactive MD

## Theoretical Methods
- Differentiable ReaxFF (TensorFlow)
- ReaxFF-MPNN hybrid potential
- Gradient-based optimization
- Neural network augmentation
- Higher-order derivatives

## Capabilities (CRITICAL)
- Differentiable ReaxFF
- ReaxFF-MPNN hybrid
- Gradient-based fitting
- Neural network augmentation
- LAMMPS integration

**Sources**: GitHub repository, Phys. Chem. Chem. Phys. 23, 19457 (2021)

## Key Strengths

### Differentiable:
- Exact gradients
- Higher-order derivatives
- TensorFlow backend
- End-to-end optimization

### Hybrid:
- ReaxFF + MPNN
- Best of both worlds
- Reactive chemistry
- ML accuracy boost

### Integration:
- LAMMPS compatible
- ReaxFF-nn for LAMMPS
- Standard workflow
- Python interface

## Inputs & Outputs
- **Input formats**: ReaxFF parameters, training data
- **Output data types**: Optimized ReaxFF parameters, hybrid models

## Interfaces & Ecosystem
- **TensorFlow**: Backend
- **LAMMPS**: MD engine
- **Python**: Core

## Performance Characteristics
- **Speed**: Standard ReaxFF + NN overhead
- **Accuracy**: Better than standard ReaxFF
- **System size**: Any (force field)
- **Automation**: Semi-automated

## Computational Cost
- **Fitting**: Hours
- **MD**: Similar to ReaxFF

## Limitations & Known Constraints
- **TensorFlow**: Required
- **Complex setup**: Multi-step
- **Limited documentation**: Academic code
- **ReaxFF basis**: Inherited limitations

## Comparison with Other Codes
- **vs JAX-ReaxFF**: I-ReaxFF adds NN, JAX-ReaxFF is optimization only
- **vs standard ReaxFF**: I-ReaxFF is differentiable + NN
- **vs pure MLIP**: I-ReaxFF has ReaxFF physics
- **Unique strength**: Differentiable ReaxFF with ReaxFF-MPNN hybrid combining reactive FF with neural networks

## Application Areas

### Reactive MD:
- Combustion
- Catalysis
- Battery electrolytes
- Polymer degradation

### Hybrid Potential:
- ReaxFF + ML accuracy
- Reactive + data-driven
- Transfer learning

## Best Practices
- Start from existing ReaxFF parameters
- Use gradient-based optimization
- Validate with MD benchmarks
- Compare with standard ReaxFF

## Community and Support
- Open source
- GitHub repository
- Academic development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/fenggo/I-ReaxFF

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Differentiable ReaxFF with ReaxFF-MPNN hybrid combining reactive FF with neural networks
