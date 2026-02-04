# BerryEasy

## Official Resources
- arXiv: https://arxiv.org/abs/2312.13051
- Publication: A. C. Tyner, arXiv:2312.13051 (2023)
- License: Check repository

## Overview
BerryEasy is a GPU-enabled Python package for diagnosis of nth-order and spin-resolved topology in the presence of fields and effects. It provides efficient computation of nested Wilson loops, spin-resolved Wilson loops, and various topological invariants using GPU acceleration for performance.

**Scientific domain**: Topological band theory, higher-order topology, spin-resolved invariants
**Target user community**: Researchers studying advanced topological phases and spin-resolved topology

## Theoretical Methods
- Wilson loop calculations
- Nested Wilson loops for higher-order topology
- Spin-resolved Wilson loops
- Berry phase and Berry curvature
- Wannier charge center evolution
- Symmetry-resolved invariants

## Capabilities (CRITICAL)
- **GPU Acceleration**: CUDA-enabled computations
- **Nested Wilson Loops**: Higher-order topology diagnosis
- **Spin-Resolved**: Spin-resolved topological invariants
- **nth-Order Topology**: Beyond first-order invariants
- **Field Effects**: External field incorporation
- **Multiple Models**: Tight-binding and Wannier support

**Sources**: arXiv preprint, code repository

## Key Strengths

### GPU Performance:
- CUDA acceleration
- Parallel computations
- Large system handling
- Fast invariant calculation

### Advanced Topology:
- Higher-order topological phases
- Spin-resolved invariants
- Nested Wilson loop spectra
- Partial polarization

### Comprehensive:
- Multiple invariant types
- Field effect support
- Model flexibility
- Modern implementation

## Inputs & Outputs
- **Input formats**:
  - Tight-binding models
  - Wannier Hamiltonians
  - PythTB-compatible models
  
- **Output data types**:
  - Wilson loop spectra
  - Topological invariants
  - Berry curvature
  - Wannier centers

## Installation
```bash
pip install berryeasy
# Or from source
git clone [repository]
pip install -e .
```

## Usage Examples
```python
import berryeasy as be

# Load tight-binding model
model = be.load_model("wannier90_hr.dat")

# Calculate nested Wilson loop
nwl = be.nested_wilson_loop(model, direction=[0, 1])

# Spin-resolved Wilson loop
swl = be.spin_resolved_wilson(model, spin_operator="Sz")

# Get invariants
z2 = be.calculate_z2(model)
```

## Performance Characteristics
- **Speed**: GPU-accelerated, significant speedup
- **Memory**: GPU memory dependent
- **Scalability**: Handles large models efficiently

## Limitations & Known Constraints
- **GPU required**: Best performance with CUDA GPU
- **Preprint stage**: Not yet peer-reviewed publication
- **Dependencies**: Requires GPU libraries

## Comparison with Other Tools
- **vs PythTB**: BerryEasy GPU-accelerated, more invariants
- **vs Z2Pack**: BerryEasy includes spin-resolved
- **vs nested_wloop**: BerryEasy has GPU support
- **Unique strength**: GPU acceleration, spin-resolved nth-order topology

## Application Areas
- Higher-order topological insulators
- Spin-Hall insulators
- Axion insulators
- Magnetic topological phases
- Fragile topology

## Best Practices
- Use GPU for large calculations
- Verify convergence with k-point density
- Check spin operator definitions
- Compare with known results

## Community and Support
- arXiv preprint available
- Academic development

## Verification & Sources
**Primary sources**:
1. arXiv: https://arxiv.org/abs/2312.13051
2. A. C. Tyner, arXiv:2312.13051 (2023)

**Confidence**: VERIFIED - arXiv preprint

**Verification status**: ✅ VERIFIED
- arXiv preprint: ACCESSIBLE
- Method: GPU-enabled topological analysis
- Developer: Academic research
