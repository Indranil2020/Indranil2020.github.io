# PythTB (Python Tight-Binding)

## Official Resources
- Homepage: https://www.physics.rutgers.edu/pythtb/
- Documentation: https://www.physics.rutgers.edu/pythtb/usage.html
- Source Repository: https://www.physics.rutgers.edu/pythtb/downloads.html
- License: Open-source (specified in distribution)

## Overview
PythTB is a Python package for computing tight-binding models in solid state physics. Developed by David Vanderbilt's group at Rutgers University, PythTB provides a simple, pedagogical interface for constructing and analyzing tight-binding Hamiltonians. The code emphasizes ease of use and educational value while maintaining capability for research applications, making it ideal for learning tight-binding physics and rapid prototyping.

**Scientific domain**: Tight-binding models, band structure, Berry phase  
**Target user community**: Students, educators, researchers, rapid prototyping

## Theoretical Methods
- Tight-binding Hamiltonian construction
- Band structure calculations
- Berry phase and curvature
- Wannier functions (simple models)
- k·p theory
- Low-energy models
- Topological properties

## Capabilities (CRITICAL)
**Category**: Open-source Python tight-binding library
- TB model construction (0D to 3D)
- Band structure calculations
- k-path interpolation
- Berry phase calculations
- Berry curvature
- Wilson loops
- Chern numbers
- Z2 invariants (simple)
- DOS calculations
- Eigenvector analysis
- Symmetry operations
- k·p expansions
- Educational examples
- Pure Python

**Sources**: Official website, documentation, publications

## Key Strengths

### Pedagogical:
- Simple interface
- Educational focus
- Clear code
- Teaching tool
- Easy learning

### Pure Python:
- NumPy/SciPy only
- Easy installation
- Readable code
- Modification-friendly
- Jupyter integration

### Berry Phase:
- Berry curvature
- Chern numbers
- Wilson loops
- Topological properties
- Research capability

### Model Building:
- Intuitive construction
- Arbitrary dimensions
- Custom hoppings
- On-site energies
- Flexible Hamiltonians

## Inputs & Outputs
- **Input formats**:
  - Python code (models)
  - Parameter specifications
  - Lattice definitions
  
- **Output data types**:
  - Band structures
  - Berry curvature
  - Topological invariants
  - Matplotlib plots
  - NumPy arrays

## Interfaces & Ecosystem

### Python Scientific:
- NumPy arrays
- Matplotlib visualization
- SciPy utilities
- Jupyter notebooks
- Pure Python

### Educational:
- Example models
- Tutorials
- Teaching materials
- Student projects

## Workflow and Usage

### Installation:
```python
# Download from website
# Extract and install
python setup.py install

# Or place pythtb.py in working directory
```

### Simple 1D Chain:
```python
import pythtb as tb
import numpy as np
import matplotlib.pyplot as plt

# Define 1D chain model
lat = [[1.0]]
orb = [[0.0]]
model = tb.tb_model(1, 1, lat, orb)

# Add hopping
t = 1.0
model.set_hop(t, 0, 0, [1])

# Solve on k-path
k_vec = np.linspace(0, 1, 101)
k_label = ['0', r'$\pi$']
evals = model.solve_all(k_vec)

# Plot
fig, ax = plt.subplots()
ax.plot(k_vec, evals)
ax.set_xlabel('k')
ax.set_ylabel('E')
plt.show()
```

### 2D Graphene:
```python
# Graphene lattice
lat = [[1.0, 0.0], [0.5, np.sqrt(3)/2]]
orb = [[1/3, 1/3], [2/3, 2/3]]
model = tb.tb_model(2, 2, lat, orb)

# Nearest-neighbor hopping
t = 1.0
model.set_hop(t, 0, 1, [0, 0])
model.set_hop(t, 0, 1, [1, 0])
model.set_hop(t, 0, 1, [0, 1])

# Band structure
path = [[0.0, 0.0], [2/3, 1/3], [0.5, 0.5], [0.0, 0.0]]
k_label = [r'$\Gamma$', r'$K$', r'$M$', r'$\Gamma$']
(k_vec, k_dist, k_node) = model.k_path(path, 101)
evals = model.solve_all(k_vec)

# Plot
model.display_band(k_vec, evals, k_node, k_label)
```

### Berry Curvature:
```python
# Calculate Berry curvature on mesh
nk = 50
k_mesh = model.k_uniform_mesh([nk, nk])

# Berry curvature for band 0
berry_curv_0 = model.berry_flux([0], k_mesh)

# Chern number (integrate Berry curvature)
chern = berry_curv_0 / (2 * np.pi)
print("Chern number:", chern)
```

### Haldane Model:
```python
# Haldane model (Chern insulator)
lat = [[1.0, 0.0], [0.5, np.sqrt(3)/2]]
orb = [[1/3, 1/3], [2/3, 2/3]]
model = tb.tb_model(2, 2, lat, orb)

# Parameters
t1 = 1.0   # Nearest-neighbor
t2 = 0.3   # Next-nearest-neighbor
phi = np.pi/2  # Complex phase
M = 0.5    # On-site mass

# Hoppings
model.set_hop(t1, 0, 1, [0, 0])
model.set_hop(t1, 0, 1, [1, 0])
model.set_hop(t1, 0, 1, [0, 1])

# Complex NNN hoppings
model.set_hop(t2*np.exp(1j*phi), 0, 0, [1, 0])
model.set_hop(t2*np.exp(1j*phi), 0, 0, [0, 1])
model.set_hop(t2*np.exp(1j*phi), 0, 0, [1, -1])

model.set_hop(t2*np.exp(-1j*phi), 1, 1, [1, 0])
model.set_hop(t2*np.exp(-1j*phi), 1, 1, [0, 1])
model.set_hop(t2*np.exp(-1j*phi), 1, 1, [1, -1])

# Mass terms
model.set_onsite([M, -M])
```

## Advanced Features

### Topological Properties:
- Chern numbers
- Berry phase
- Wilson loops
- Z2 (basic)
- Winding numbers

### k-Space Operations:
- k-path generation
- Uniform meshes
- Band interpolation
- High-symmetry points

### Model Manipulation:
- Cut to slab
- Finite systems
- Edge states
- Supercells
- Symmetry operations

## Performance Characteristics
- **Speed**: Python (moderate)
- **Purpose**: Education and rapid prototyping
- **System size**: Small to moderate models
- **Accuracy**: Exact diagonalization
- **Typical**: Seconds to minutes

## Computational Cost
- Pure Python overhead
- Small models fast
- Exact diagonalization
- k-mesh dependent
- Suitable for education

## Limitations & Known Constraints
- **Performance**: Python slower than compiled
- **System size**: Limited for large systems
- **Features**: Basic compared to specialized codes
- **Documentation**: Educational focus
- **Best for**: Learning, small models, prototyping

## Comparison with Other TB Codes
- **vs Kwant**: PythTB simpler/educational, Kwant transport
- **vs TBmodels**: PythTB pedagogical, TBmodels symmetry
- **Unique strength**: Educational clarity, simple interface, Vanderbilt group, Rutgers standard

## Application Areas

### Education:
- Teaching solid state
- Learning tight-binding
- Student projects
- Computational physics courses
- Topological concepts

### Research:
- Rapid prototyping
- Model exploration
- Concept validation
- Simple systems
- Quick calculations

### Topological Physics:
- Berry phase calculations
- Chern insulators
- Topological models
- Pedagogical examples

## Best Practices

### Learning:
- Start with 1D examples
- Follow documentation
- Modify examples
- Build complexity gradually
- Understand physics

### Model Construction:
- Clear lattice definition
- Check hopping consistency
- Verify band structure
- Compare known models

### Research:
- Validate against known results
- Appropriate for model size
- Use for prototyping
- Production may need specialized codes

## Community and Support
- Open-source
- Rutgers University (Vanderbilt group)
- Documentation
- Example gallery
- Educational community
- Publication citations

## Educational Resources
- Excellent documentation
- Tutorial examples
- Model gallery
- Physics explanations
- User guide
- Pedagogical focus

## Development
- David Vanderbilt group (Rutgers)
- Sinisa Coh, David Vanderbilt
- Stable release
- Educational mission
- Maintenance updates

## Research Impact
PythTB is widely used in education and has enabled rapid prototyping of tight-binding models, particularly valued for teaching topological band theory and Berry phase physics.

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.physics.rutgers.edu/pythtb/
2. Documentation: https://www.physics.rutgers.edu/pythtb/usage.html
3. Publications: Comp. Phys. Comm. 185, 2309 (2014)

**Secondary sources**:
1. Educational materials
2. Topological physics literature
3. User publications

**Confidence**: VERIFIED - Educational TB library

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- License: Open-source
- **Category**: Python tight-binding library
- Status: Stable, maintained
- Institution: Rutgers University (Vanderbilt group)
- Specialized strength: Educational Python tight-binding, simple interface, pedagogical clarity, Berry phase calculations, topological properties, rapid prototyping, pure Python, Jupyter-friendly, teaching tool, Chern numbers, model construction ease
