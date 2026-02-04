# pysktb

## Official Resources
- **GitHub**: https://github.com/santoshkumarradha/pysktb
- **Documentation**: Available in repository
- **License**: MIT License

## Overview
pysktb is a Python package for Slater-Koster tight-binding calculations with a focus on topological materials analysis. It provides tools for constructing tight-binding Hamiltonians using the two-center Slater-Koster approximation and calculating topological invariants like Berry phases and Chern numbers.

**Scientific domain**: Tight-binding models, topological materials, electronic structure
**Target user community**: Researchers studying topological materials, 2D systems, and model Hamiltonians

## Theoretical Background
pysktb implements:
- Slater-Koster two-center approximation: H_ij = Σ V_llm(r) Y_lm(r̂)
- Spin-orbit coupling: H_SOC = λ L·S
- Berry phase: γ = i∮⟨u_k|∇_k|u_k⟩dk
- Chern number: C = (1/2π)∫ F dk
- Z2 topological invariants

## Capabilities (CRITICAL)
- **Slater-Koster TB**: Two-center approximation for s, p, d orbitals
- **Topological Analysis**: Berry phase, Chern numbers, Z2 invariants
- **Band Structure**: Electronic band calculations along k-paths
- **Spin-Orbit Coupling**: Full SOC implementation
- **Edge States**: Surface/edge state calculations
- **Wannier Centers**: Hybrid Wannier function analysis

## Key Strengths

### Slater-Koster Framework:
- Standard two-center parameters
- s, p, d orbital support
- Customizable hopping parameters
- Distance-dependent interactions

### Topological Analysis:
- Berry phase calculation
- Chern number computation
- Z2 invariant determination
- Wilson loop analysis

### Spin-Orbit Coupling:
- Atomic SOC implementation
- Rashba and Dresselhaus terms
- Spin texture visualization

## Inputs & Outputs
- **Input formats**:
  - Python dictionaries for parameters
  - Structure definitions
  - Hopping parameters
  
- **Output data types**:
  - Band structures
  - Topological invariants
  - Berry curvature
  - Edge state spectra

## Installation
```bash
pip install pysktb
```

## Usage Examples
```python
from pysktb import Lattice, TightBinding

# Define lattice
lattice = Lattice([[1,0],[0,1]])

# Create tight-binding model
tb = TightBinding(lattice)
tb.add_orbital([0,0], 's')
tb.add_hopping(1, 0, 1, -1.0)  # Nearest neighbor hopping

# Calculate band structure
kpath = [[0,0], [0.5,0], [0.5,0.5], [0,0]]
bands = tb.solve_along_path(kpath, 100)

# Calculate Chern number
chern = tb.chern_number(band_index=0)
```

## Performance Characteristics
- **Speed**: Fast for model Hamiltonians
- **Memory**: Efficient sparse matrices
- **Scalability**: Suitable for large unit cells

## Limitations & Known Constraints
- **Model-based**: Requires parameter input
- **Two-center**: Limited to Slater-Koster approximation
- **Parameterization**: User must provide hopping parameters

## Comparison with Other Tools
- **vs PythTB**: Similar capabilities, different API
- **vs sisl**: pysktb focused on topological analysis
- **vs Kwant**: pysktb simpler for bulk calculations
- **Unique strength**: Topological invariant calculations

## Application Areas
- Topological insulators
- Weyl/Dirac semimetals
- 2D materials (graphene, TMDs)
- Heterostructures
- Quantum spin Hall systems
- Topological crystalline insulators

## Best Practices
- Validate parameters against DFT
- Check topological invariant convergence
- Use appropriate k-point sampling
- Verify edge state localization

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/santoshkumarradha/pysktb

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Developer: Santosh Kumar Radha
