# Pybinding

## Official Resources
- Homepage: https://docs.pybinding.site/
- Documentation: https://docs.pybinding.site/en/stable/
- Source Repository: https://github.com/dean0x7d/pybinding
- License: BSD 3-Clause License

## Overview
Pybinding is a Python package for tight-binding calculations with a focus on large-scale systems and efficient computation. Developed by Dean Moldovan, pybinding provides a high-level interface for constructing tight-binding models with C++ performance under the hood. The code emphasizes computational efficiency, enabling large-scale calculations while maintaining Python's ease of use for model construction and analysis.

**Scientific domain**: Tight-binding, large-scale systems, computational efficiency  
**Target user community**: Python users, large systems, performance-focused calculations

## Theoretical Methods
- Tight-binding Hamiltonian
- Sparse matrix techniques
- Band structure calculations
- LDOS and spectral functions
- Green's function methods
- Kernel polynomial method (KPM)
- Large-scale eigensolvers

## Capabilities (CRITICAL)
**Category**: Open-source Python TB library
- TB model construction
- Efficient C++ core
- Large-scale systems (10⁶+ sites)
- Band structure calculations
- Density of states (KPM)
- Spectral functions
- LDOS calculations
- Arbitrary geometries
- Periodic/finite systems
- Magnetic fields
- Disorder
- Parallel computation
- Visualization tools
- Production quality

**Sources**: Official documentation, GitHub

## Key Strengths

### Performance:
- C++ core implementation
- Sparse matrix optimization
- Large-scale capable
- KPM efficiency
- Production speed

### Python Interface:
- High-level API
- Easy model building
- NumPy integration
- Visualization
- User-friendly

### KPM Implementation:
- Efficient DOS
- Large systems
- Spectral functions
- O(N) scaling
- Production quality

### Flexibility:
- Arbitrary lattices
- Complex geometries
- Custom modifiers
- Research capable

## Inputs & Outputs
- **Input formats**:
  - Python model definition
  - Lattice specifications
  - Parameters
  
- **Output data types**:
  - Band structures
  - DOS
  - LDOS
  - Spectral functions
  - NumPy arrays
  - Visualization plots

## Interfaces & Ecosystem

### Python Scientific:
- NumPy/SciPy
- Matplotlib
- Jupyter notebooks
- Modern Python

### Visualization:
- Built-in plotting
- System visualization
- Band structure plots
- LDOS maps

## Workflow and Usage

### Installation:
```bash
# Via conda (recommended)
conda install -c conda-forge pybinding

# Via pip
pip install pybinding
```

### Simple Graphene:
```python
import pybinding as pb
import numpy as np
import matplotlib.pyplot as plt

# Define lattice
lattice = pb.Lattice(
    a1=[1, 0],
    a2=[0.5, 0.5 * np.sqrt(3)]
)

# Add sublattices
lattice.add_sublattices(
    ('A', [0, 0]),
    ('B', [0, 1/np.sqrt(3)])
)

# Add hoppings
lattice.add_hoppings(
    ([0, 0], 'A', 'B', -1),
    ([1, 0], 'A', 'B', -1),
    ([0, 1], 'A', 'B', -1)
)

# Build model
model = pb.Model(lattice)

# Compute band structure
solver = pb.solver.lapack(model)
bands = solver.calc_bands(-np.pi, np.pi)
bands.plot()
plt.show()
```

### Finite System:
```python
# Graphene nanoribbon
model = pb.Model(
    lattice,
    pb.rectangle(10, 10)  # 10x10 nm rectangle
)

# Solve
solver = pb.solver.arpack(model, k=20)  # 20 eigenvalues
eigenvalues = solver.eigenvalues
eigenvectors = solver.eigenvectors
```

### Density of States (KPM):
```python
# Large system DOS using KPM
model = pb.Model(
    lattice,
    pb.rectangle(100, 100)  # Large system
)

# KPM calculation
kpm = pb.kpm(model)
dos = kpm.calc_dos(
    energy=np.linspace(-3, 3, 500),
    broadening=0.05,
    num_moments=1000
)

dos.plot()
plt.show()
```

### LDOS Calculation:
```python
# Local density of states
ldos = kpm.calc_ldos(
    energy=np.linspace(-1, 1, 100),
    broadening=0.05,
    position=[0, 0]  # Position in system
)

ldos.plot()
```

### Magnetic Field:
```python
# Magnetic field via Peierls substitution
@pb.hopping_modifier
def constant_magnetic_field(energy, hop, B=0):
    x, y = hop.r[:, 0], hop.r[:, 1]
    phase = 2 * np.pi * B * x * y / 2
    return energy * np.exp(1j * phase)

model = pb.Model(
    lattice,
    pb.rectangle(20, 20),
    constant_magnetic_field(B=10)
)
```

### Disorder:
```python
# Random on-site disorder
@pb.onsite_energy_modifier
def random_potential(energy, sub_id):
    return energy + np.random.uniform(-0.5, 0.5, energy.shape)

model = pb.Model(
    lattice,
    pb.rectangle(10, 10),
    random_potential
)
```

## Advanced Features

### KPM Method:
- Kernel polynomial method
- Efficient large systems
- DOS and spectral functions
- LDOS calculations
- Production quality

### Modifiers:
- Onsite energy modifiers
- Hopping modifiers
- Magnetic fields
- Disorder
- Custom potentials

### Solvers:
- LAPACK (small systems)
- ARPACK (sparse eigenvalues)
- FEAST (energy intervals)
- KPM (large DOS)
- Appropriate selection

### Parallel:
- OpenMP threading
- Multi-core support
- Large-scale performance
- Efficient computation

## Performance Characteristics
- **Speed**: C++ core, very fast
- **System size**: 10⁶+ sites (KPM)
- **Purpose**: Large-scale TB
- **Scalability**: Excellent
- **Typical**: Seconds to minutes

## Computational Cost
- C++ efficiency
- Sparse matrix optimized
- KPM O(N) scaling
- Memory efficient
- Production capable

## Limitations & Known Constraints
- **Lattice-based**: Discrete systems
- **Python overhead**: Interface layer
- **Documentation**: Good but evolving
- **Maintenance**: Single developer primarily
- **Transport**: Not transport-focused (see Kwant)

## Comparison with Other TB Codes
- **vs Kwant**: Pybinding performance/DOS, Kwant transport
- **vs pythtb**: Pybinding large systems, pythtb pedagogical
- **Unique strength**: C++ performance with Python ease, KPM efficiency, large-scale capable, modern design

## Application Areas

### Large-Scale Systems:
- Graphene devices
- Large nanostructures
- Disordered systems
- Statistical sampling
- Realistic sizes

### Spectral Properties:
- DOS calculations
- LDOS maps
- Spectral functions
- Large system spectra
- Disorder averaging

### 2D Materials:
- Graphene
- TMDs
- Heterostructures
- Moiré patterns
- Device modeling

## Best Practices

### Model Construction:
- Define lattice carefully
- Use modifiers efficiently
- Test on small systems
- Validate Hamiltonians

### Large Systems:
- Use KPM for DOS
- Sparse eigensolvers
- Memory awareness
- Parallel execution

### Performance:
- Appropriate solver choice
- KPM parameters
- Threading
- Profiling

## Community and Support
- Open-source (BSD 3-Clause)
- GitHub repository
- Documentation
- Active development
- Issue tracking
- User community

## Educational Resources
- Comprehensive documentation
- Tutorial examples
- Example gallery
- API reference
- Physics explanations

## Development
- Dean Moldovan (lead)
- Active development
- Regular updates
- Community contributions
- Performance focus

## Verification & Sources
**Primary sources**:
1. Homepage: https://docs.pybinding.site/
2. GitHub: https://github.com/dean0x7d/pybinding
3. Documentation

**Secondary sources**:
1. User publications
2. TB literature

**Confidence**: VERIFIED - Performance TB library

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: BSD 3-Clause (open-source)
- **Category**: Python TB library
- Status: Actively developed
- Specialized strength: High-performance tight-binding with Python interface, C++ core, large-scale systems (10⁶+ sites), kernel polynomial method (KPM), efficient DOS/LDOS calculations, sparse matrix optimization, modern design, production quality, computational efficiency
