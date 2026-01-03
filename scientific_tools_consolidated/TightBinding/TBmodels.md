# TBmodels

## Official Resources
- Homepage: https://tbmodels.greschd.ch/
- Documentation: https://tbmodels.greschd.ch/en/latest/
- Source Repository: https://github.com/Z2PackDev/TBmodels
- License: GNU General Public License v3.0

## Overview
TBmodels is a Python library for evaluating tight-binding models with an emphasis on symmetry and topological properties. Developed by Dominik Gresch as part of the Z2Pack ecosystem, TBmodels provides tools for constructing, manipulating, and analyzing tight-binding Hamiltonians with automatic symmetry operations. The code integrates seamlessly with Wannier90 and Z2Pack, enabling workflows from ab-initio calculations to topological characterization.

**Scientific domain**: Tight-binding models, symmetry operations, topological physics  
**Target user community**: Python users, symmetry-aware calculations, Z2Pack users

## Theoretical Methods
- Tight-binding Hamiltonian evaluation
- Symmetry operations
- k-space interpolation
- Wannier function transformations
- Band structure calculations
- Symmetry-adapted basis

## Capabilities (CRITICAL)
**Category**: Open-source Python TB library
- TB model construction and manipulation
- Symmetry operations (automatic)
- Wannier90 integration
- k-space evaluation
- Band structure calculations
- Model transformations
- Sparse representations
- Real-space to k-space
- Z2Pack integration
- HDF5 I/O
- Production quality

**Sources**: Official documentation, GitHub

## Key Strengths

### Symmetry-Aware:
- Automatic symmetry operations
- Symmetrization
- Symmetry-adapted basis
- Efficient implementation
- Research quality

### Z2Pack Ecosystem:
- Native integration
- Topological workflows
- Seamless pipeline
- Wannier90 connection
- Modern tools

### Python Clean Design:
- Object-oriented
- Well-documented
- Type hints
- Modern Python
- Extensible

### Wannier90 Interface:
- Direct hr.dat reading
- Transformation tools
- Standard workflow
- Production integration

## Inputs & Outputs
- **Input formats**:
  - Wannier90 hr.dat
  - HDF5 files
  - Python objects
  - Manual construction
  
- **Output data types**:
  - TB model objects
  - Band structures
  - HDF5 storage
  - NumPy arrays
  - Symmetrized models

## Interfaces & Ecosystem

### Wannier90:
- hr.dat import
- Model conversion
- Transformation tools

### Z2Pack:
- Native integration
- Topological calculations
- Workflow compatibility

### Python Scientific:
- NumPy/SciPy
- HDF5
- Matplotlib
- Modern ecosystem

## Workflow and Usage

### Installation:
```bash
pip install tbmodels
```

### Load Wannier90 Model:
```python
import tbmodels as tb

# Load from Wannier90
model = tb.Model.from_wannier_files(
    hr_file='wannier90_hr.dat',
    wsvec_file='wannier90_wsvec.dat'
)

# Or from hr.dat only
model = tb.Model.from_hr_file('wannier90_hr.dat')
```

### Evaluate at k-points:
```python
import numpy as np

# Single k-point
k = [0.0, 0.0, 0.0]
hamiltonian = model.hamilton(k)
eigenvalues = np.linalg.eigvalsh(hamiltonian)

# k-path
k_points = np.array([[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0]])
bands = [np.linalg.eigvalsh(model.hamilton(k)) for k in k_points]
```

### Symmetry Operations:
```python
from tbmodels import symmetry

# Define symmetry operation
rotation = symmetry.rotation([0, 0, 1], 3)  # C3 rotation

# Apply to model
model_sym = symmetry.symmetrize(
    model,
    symmetries=[rotation]
)
```

### Band Structure:
```python
# Generate k-path
k_points = tb.kpoints.path(
    ['Gamma', 'X', 'M', 'Gamma'],
    lattice=model.uc
)

# Calculate bands
eigenvalues = np.array([
    np.linalg.eigvalsh(model.hamilton(k))
    for k in k_points
])

# Plot
import matplotlib.pyplot as plt
plt.plot(eigenvalues)
plt.ylabel('Energy (eV)')
plt.show()
```

### Z2Pack Integration:
```python
import z2pack

# Create Z2Pack system from TBmodels
system = z2pack.tb.System(model)

# Run topological calculation
result = z2pack.surface.run(
    system=system,
    surface=lambda s, t: [s/2, t, 0]
)
```

### Save/Load:
```python
# Save to HDF5
model.to_hdf5_file('model.hdf5')

# Load
model_loaded = tb.Model.from_hdf5_file('model.hdf5')
```

## Advanced Features

### Model Manipulation:
- Slice orbitals
- Add on-site terms
- Modify hoppings
- Transform coordinates
- Real-space construction

### Sparse Representations:
- Memory efficient
- Large systems
- Fast evaluation
- Production capable

### Symmetrization:
- Automatic symmetry
- Custom operations
- Basis transformation
- High-symmetry enforcement

## Performance Characteristics
- **Speed**: Efficient Python
- **Purpose**: TB model manipulation
- **System size**: Any (post-Wannier90)
- **Memory**: Efficient sparse
- **Typical**: Interactive to production

## Computational Cost
- Post-processing tool
- Fast evaluation
- k-point dependent
- Memory efficient
- Production suitable

## Limitations & Known Constraints
- **Python overhead**: Not fastest
- **Requires input**: Wannier90 or manual
- **Scope**: TB model tool, not DFT
- **Documentation**: Good but Z2Pack-focused

## Comparison with Other TB Libraries
- **vs pythtb**: TBmodels symmetry-focused, pythtb pedagogical
- **vs Kwant**: TBmodels general TB, Kwant transport
- **Unique strength**: Symmetry operations, Z2Pack integration, Wannier90 interface, modern Python design

## Application Areas

### Topological Calculations:
- Z2Pack workflows
- Symmetry-aware topology
- Wannier-based topology
- High-throughput

### Symmetry Analysis:
- Symmetrization
- Symmetry operations
- Irreducible representations
- Crystal symmetries

### Model Manipulation:
- Wannier90 post-processing
- Model transformations
- Band structure analysis
- Research workflows

## Best Practices

### Wannier90 Input:
- Quality MLWFs
- Proper wsvec file
- Validated models
- Standard workflow

### Symmetry:
- Define operations correctly
- Check symmetrization results
- Use appropriate group
- Validate bands

### Integration:
- Z2Pack for topology
- Combine with WannierTools
- Python ecosystem
- Modular usage

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Documentation
- Z2Pack ecosystem
- Active development

## Educational Resources
- API documentation
- Example gallery
- Z2Pack tutorials
- Symmetry guides

## Development
- Dominik Gresch (lead)
- Z2Pack developers
- Active maintenance
- Regular updates
- Community contributions

## Verification & Sources
**Primary sources**:
1. Homepage: https://tbmodels.greschd.ch/
2. GitHub: https://github.com/Z2PackDev/TBmodels
3. Z2Pack documentation

**Secondary sources**:
1. Z2Pack publications
2. User applications

**Confidence**: VERIFIED - Python TB library

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: GPL v3 (open-source)
- **Category**: Python tight-binding library
- Status: Actively maintained
- Specialized strength: Symmetry-aware tight-binding models, automatic symmetry operations, Wannier90 integration, Z2Pack ecosystem, model manipulation, HDF5 I/O, modern Python design, sparse representations, topological workflows
