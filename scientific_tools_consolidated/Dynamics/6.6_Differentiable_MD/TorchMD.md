# TorchMD

## Official Resources
- Homepage: https://github.com/torchmd/torchmd
- Documentation: https://torchmd.readthedocs.io/
- Source Repository: https://github.com/torchmd/torchmd
- License: MIT

## Overview
TorchMD is an end-to-end molecular dynamics engine using PyTorch. It enables differentiable simulations, seamless integration of neural network potentials, and end-to-end training of ML models through molecular dynamics trajectories.

**Scientific domain**: Differentiable MD, PyTorch-based simulations, ML potentials  
**Target user community**: Researchers integrating ML with molecular dynamics

## Theoretical Methods
- Differentiable molecular dynamics
- Neural network potentials
- PyTorch automatic differentiation
- Classical force fields
- End-to-end learning

## Capabilities (CRITICAL)
- End-to-end differentiable MD
- PyTorch backend
- Neural network potential integration
- GPU acceleration
- Force field support
- TorchMD-NET integration

## Key Strengths

### PyTorch Integration:
- Native PyTorch
- Easy ML integration
- Automatic differentiation
- GPU acceleration

### ML Potentials:
- TorchMD-NET
- Custom potentials
- End-to-end training

## Inputs & Outputs
- **Input formats**:
  - PDB structures
  - PyTorch tensors
  - Parameter files
  
- **Output data types**:
  - Trajectories
  - Gradients
  - Energies/forces

## Interfaces & Ecosystem
- **PyTorch**: Backend
- **TorchMD-NET**: Neural potentials
- **OpenMM**: Compatibility

## Advanced Features
- **Differentiable**: Full autodiff
- **TorchMD-NET**: Equivariant NNs
- **Custom potentials**: Easy to add
- **GPU acceleration**: CUDA support

## Performance Characteristics
- Good GPU performance
- PyTorch efficiency
- Autodiff overhead
- Good for ML integration

## Computational Cost
- GPU provides speedup
- Autodiff adds cost
- ML potentials efficient
- Overall: Good for ML workflows

## Best Practices
- Use GPU acceleration
- Validate potentials
- Check energy conservation
- Use TorchMD-NET for accuracy

## Limitations & Known Constraints
- PyTorch dependency
- Less traditional features
- Active development
- Performance vs specialized codes

## Application Areas
- ML potential development
- Drug discovery
- Coarse-grained modeling
- End-to-end learning
- Method development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/torchmd/torchmd
2. S. Doerr et al., J. Chem. Theory Comput. 17, 2355 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
