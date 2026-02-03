# JAX-MD

## Official Resources
- Homepage: https://github.com/jax-md/jax-md
- Documentation: https://jax-md.readthedocs.io/
- Source Repository: https://github.com/jax-md/jax-md
- License: Apache-2.0

## Overview
JAX-MD is a framework for differentiable physics simulations with a focus on molecular dynamics. Built on JAX, it enables end-to-end differentiable simulations, automatic differentiation through dynamics, and seamless integration with machine learning.

**Scientific domain**: Differentiable MD, machine learning, GPU-accelerated simulations  
**Target user community**: Researchers combining MD with machine learning

## Theoretical Methods
- Differentiable molecular dynamics
- Automatic differentiation
- Neural network potentials
- End-to-end learning
- Various integrators (NVE, NVT, NPT)

## Capabilities (CRITICAL)
- End-to-end differentiable MD
- JAX automatic differentiation
- GPU/TPU acceleration
- Neural network integration
- Multiple thermostats/barostats
- Custom potentials

## Key Strengths

### Differentiability:
- Full autodiff support
- Gradient through dynamics
- End-to-end training
- Loss function flexibility

### JAX Backend:
- GPU/TPU acceleration
- JIT compilation
- Vectorization (vmap)
- Functional programming

## Inputs & Outputs
- **Input formats**:
  - JAX arrays
  - Python configuration
  
- **Output data types**:
  - Trajectories
  - Gradients
  - Observables

## Interfaces & Ecosystem
- **JAX**: Backend
- **Flax/Haiku**: Neural networks
- **Optax**: Optimization

## Advanced Features
- **Autodiff**: Through entire simulation
- **Neural potentials**: Easy integration
- **Custom forces**: Differentiable
- **Thermostats**: Nose-Hoover, Langevin
- **Neighbor lists**: Efficient

## Performance Characteristics
- Excellent GPU performance
- JIT compilation
- Vectorization support
- Good scaling

## Computational Cost
- GPU provides major speedup
- Autodiff adds overhead
- JIT compilation helps
- Overall: Efficient on GPU

## Best Practices
- Use JIT compilation
- Leverage vmap for batching
- Validate gradients
- Use appropriate precision

## Limitations & Known Constraints
- JAX ecosystem required
- Different programming paradigm
- Less traditional MD features
- Active development

## Application Areas
- ML potential development
- Inverse design
- Parameter optimization
- Method development
- Differentiable physics

## Comparison with Other Codes
- **vs TorchMD**: JAX-MD uses JAX, TorchMD uses PyTorch
- **vs OpenMM**: JAX-MD fully differentiable, OpenMM traditional MD
- **vs LAMMPS**: JAX-MD for ML research, LAMMPS for production
- **Unique strength**: Full autodiff through dynamics, JAX ecosystem (vmap, pmap, JIT), GPU/TPU

## Community and Support
- Google Research development
- GitHub issues
- Growing community
- JAX ecosystem

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/jax-md/jax-md
2. S. Schoenholz & E. Cubuk, NeurIPS 2020

**Secondary sources**:
1. JAX-MD tutorials
2. JAX ecosystem documentation
3. Published ML physics applications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, Apache-2.0)
- Academic citations: >400
- Active development: Google Research
- Growing adoption in ML physics
