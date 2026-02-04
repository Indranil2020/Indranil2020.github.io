# topo_2bands

## Official Resources
- GitHub: https://github.com/jameclear/topo_2bands
- License: Check repository

## Overview
topo_2bands is a program for calculating typical topological invariants based on minimum two-band tight-binding models. It computes Berry phase, winding number, Chern number, and Z2 invariant, as well as edge states and Fermi arc states in 2D and 3D gapless systems.

**Scientific domain**: Topological band theory, tight-binding models, topological invariants
**Target user community**: Researchers and students learning topological band theory with model Hamiltonians

## Theoretical Methods
- Berry phase calculation
- Winding number computation
- Chern number integration
- Z2 invariant determination
- Edge state calculation
- Fermi arc visualization

## Capabilities (CRITICAL)
- **Berry Phase**: Single band Berry phase
- **Winding Number**: 1D topological invariant
- **Chern Number**: 2D topological invariant
- **Z2 Invariant**: Time-reversal protected topology
- **Edge States**: Boundary mode calculation
- **Fermi Arcs**: 3D surface state arcs

**Sources**: GitHub repository

## Key Strengths

### Educational Value:
- Simple two-band models
- Clear implementations
- Learning-focused
- Example gallery

### Multiple Invariants:
- Comprehensive coverage
- Standard definitions
- Validated results
- Cross-comparisons

### Visualization:
- Edge state plotting
- Fermi arc display
- Band structure visualization
- Invariant evolution

## Inputs & Outputs
- **Input formats**:
  - Model parameters
  - Hamiltonian definition
  - k-space discretization
  
- **Output data types**:
  - Topological invariants
  - Edge state spectra
  - Fermi arc plots
  - Berry phase evolution

## Installation
```bash
git clone https://github.com/jameclear/topo_2bands.git
cd topo_2bands
# Follow installation instructions
```

## Usage Examples
```python
# Example: Calculate Chern number for two-band model
from topo_2bands import chern_number, edge_states

# Define model parameters
params = {
    'm': 1.0,    # Mass term
    't': 1.0,    # Hopping
    'lambda_so': 0.5  # SOC
}

# Calculate Chern number
C = chern_number(params, nk=100)
print(f"Chern number: {C}")

# Calculate edge states
E_edge, k_edge = edge_states(params, Ly=50)
```

## Performance Characteristics
- **Speed**: Fast for two-band models
- **Accuracy**: Numerical precision controlled
- **Simplicity**: Easy to understand and modify

## Limitations & Known Constraints
- **Two-band only**: Limited to two-band models
- **Model-based**: Not for ab initio calculations
- **Educational focus**: Not production code

## Comparison with Other Tools
- **vs PythTB**: topo_2bands simpler, specialized
- **vs Z2Pack**: topo_2bands educational focus
- **vs pyqula**: topo_2bands more focused scope
- **Unique strength**: Clear educational implementation

## Application Areas
- Topological band theory education
- Model Hamiltonian studies
- Invariant computation learning
- Research prototyping
- Method development

## Best Practices
- Start with known models (SSH, Haldane)
- Verify against analytical results
- Check convergence with k-points
- Use for understanding before applying

## Community and Support
- GitHub repository
- Educational use
- Example-driven

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/jameclear/topo_2bands

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Source code: OPEN
- Method: Two-band topological invariants
