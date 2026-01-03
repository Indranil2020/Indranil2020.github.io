# Kwant

## Official Resources
- Homepage: https://kwant-project.org/
- Documentation: https://kwant-project.org/doc/
- Source Repository: https://gitlab.kwant-project.org/kwant/kwant
- License: BSD 2-Clause License

## Overview
Kwant is a powerful Python package for numerical quantum transport calculations. Developed as a community project with contributions from multiple institutions, Kwant enables the simulation of tight-binding systems with an emphasis on mesoscopic physics, quantum transport, and topological phenomena. The code provides an intuitive interface for constructing complex geometries, calculating conductance, and analyzing scattering properties in low-dimensional systems.

**Scientific domain**: Quantum transport, mesoscopic physics, tight-binding  
**Target user community**: Mesoscopic physics, quantum transport, nanostructures, topological transport

## Theoretical Methods
- Tight-binding Hamiltonian construction
- Landauer-Büttiker formalism
- Scattering matrix calculations
- Green's function methods
- Recursive Green's functions
- Lead self-energies
- Wave function calculations
- Local density of states

## Capabilities (CRITICAL)
**Category**: Open-source quantum transport library
- TB system construction (flexible geometry)
- Conductance calculations
- Scattering matrix (S-matrix)
- Wave functions
- Local density of states (LDOS)
- Arbitrary system shapes
- Multiple leads
- Magnetic fields
- Disorder and impurities
- Superconductivity (BdG)
- Topological transport
- Parallel execution
- Python implementation
- Production quality

**Sources**: Official website, documentation, publications

## Key Strengths

### Quantum Transport:
- Landauer formalism
- Multi-terminal systems
- Conductance calculations
- Scattering properties
- Industry standard

### Geometry Flexibility:
- Arbitrary shapes
- Lattice-based construction
- Builder pattern
- Complex geometries
- Intuitive interface

### Physics Coverage:
- Normal transport
- Magnetic fields
- Superconductivity (BdG)
- Topological systems
- Disorder effects

### Python Design:
- User-friendly API
- NumPy/SciPy integration
- Visualization tools
- Extensible
- Well-documented

## Inputs & Outputs
- **Input formats**:
  - Python code (system definition)
  - Lattice specifications
  - Hamiltonian parameters
  
- **Output data types**:
  - Conductance
  - S-matrix elements
  - Wave functions
  - LDOS
  - Band structures
  - NumPy arrays

## Interfaces & Ecosystem

### Python Scientific:
- NumPy arrays
- SciPy sparse matrices
- Matplotlib visualization
- Jupyter notebooks
- Modern Python

### Visualization:
- Built-in plotting
- 3D system visualization
- Wave function plots
- LDOS maps

## Workflow and Usage

### Installation:
```bash
# Via conda
conda install -c conda-forge kwant

# Via pip
pip install kwant
```

### Simple 1D Wire:
```python
import kwant
import numpy as np

# Define lattice
lat = kwant.lattice.chain(a=1)

# Build system
syst = kwant.Builder()

# Scattering region (10 sites)
for i in range(10):
    syst[lat(i)] = 0.5  # On-site energy

# Hoppings
for i in range(9):
    syst[lat(i), lat(i+1)] = -1.0

# Attach leads
lead = kwant.Builder(kwant.TranslationalSymmetry((-1,)))
lead[lat(0)] = 0.5
lead[lat(0), lat(1)] = -1.0

syst.attach_lead(lead)
syst.attach_lead(lead.reversed())

# Finalize
fsyst = syst.finalized()

# Calculate conductance
energies = np.linspace(-2, 2, 100)
conductances = [kwant.smatrix(fsyst, E).transmission(1, 0) 
                for E in energies]

# Plot
import matplotlib.pyplot as plt
plt.plot(energies, conductances)
plt.xlabel('Energy')
plt.ylabel('Conductance (e²/h)')
plt.show()
```

### 2D Quantum Hall:
```python
# Square lattice
lat = kwant.lattice.square(a=1)

# Build Hall bar
W, L = 10, 30
syst = kwant.Builder()

def shape(pos):
    x, y = pos
    return 0 <= x < L and 0 <= y < W

# Magnetic field via Peierls substitution
def hopping(site1, site2, B):
    x1, y1 = site1.pos
    x2, y2 = site2.pos
    # Landau gauge: A = (0, Bx, 0)
    phase = 2 * np.pi * B * (x1 + x2) * (y2 - y1) / 2
    return -1 * np.exp(1j * phase)

# Add sites and hoppings
syst[lat.shape(shape, (0, 0))] = 0
syst[lat.neighbors()] = hopping

# Leads
lead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))
lead[lat.shape(lambda pos: 0 <= pos[1] < W, (0, 0))] = 0
lead[lat.neighbors()] = hopping

syst.attach_lead(lead)
syst.attach_lead(lead.reversed())

# Finalize and calculate
fsyst = syst.finalized()
smat = kwant.smatrix(fsyst, energy=0.3, args=[0.01])
print("Hall conductance:", smat.transmission(1, 0))
```

### Wave Functions:
```python
# Calculate wave function
wf = kwant.wave_function(fsyst, energy=0.5)
psi = wf(0)  # Wave function from lead 0

# Plot
kwant.plotter.map(fsyst, np.abs(psi[0])**2)
```

### Graphene Nanoribbon:
```python
# Graphene lattice
graphene = kwant.lattice.general(
    [(1, 0), (0.5, np.sqrt(3)/2)],
    [(0, 0), (0, 1/np.sqrt(3))]
)
a, b = graphene.sublattices

# Build ribbon
syst = kwant.Builder()

def ribbon(pos):
    x, y = pos
    return 0 <= x < 10 and -1 <= y < 1

syst[graphene.shape(ribbon, (0, 0))] = 0
syst[graphene.neighbors()] = -1

# Add leads
lead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))
lead[graphene.shape(ribbon, (0, 0))] = 0
lead[graphene.neighbors()] = -1

syst.attach_lead(lead)
syst.attach_lead(lead.reversed())
```

## Advanced Features

### Superconductivity:
- Bogoliubov-de Gennes (BdG)
- Particle-hole space
- Andreev reflection
- Josephson junctions
- Majorana modes

### Magnetic Fields:
- Peierls substitution
- Orbital effects
- Landau levels
- Quantum Hall physics
- Aharonov-Bohm

### Disorder:
- Random on-site potentials
- Anderson localization
- Ensemble averaging
- Statistical analysis

### Sparse Matrix:
- Efficient large systems
- SciPy sparse integration
- Memory optimization
- Production capable

## Performance Characteristics
- **Speed**: Efficient Python/C
- **System size**: 10⁴-10⁶ sites
- **Purpose**: Mesoscopic transport
- **Scalability**: Good parallelization
- **Typical**: Minutes to hours

## Computational Cost
- System-size dependent
- Lead calculation important
- Sparse matrix efficient
- Energy-point parallelizable
- HPC capable

## Limitations & Known Constraints
- **Tight-binding**: Discrete lattice only
- **Leads**: Semi-infinite approximation
- **Interactions**: Mean-field
- **Temperature**: Zero-temperature formalism
- **Learning curve**: Moderate

## Comparison with Other TB Codes
- **vs pythtb**: Kwant transport focus, pythtb band structure
- **vs Pybinding**: Similar scope, Kwant more established
- **Unique strength**: Quantum transport standard, Landauer formalism, flexible geometry, production quality

## Application Areas

### Mesoscopic Physics:
- Quantum point contacts
- Quantum dots
- Nanowires
- 2D electron gas
- Hall bars

### Topological Transport:
- Edge state transport
- Quantum spin Hall
- Majorana wires
- Topological insulators
- Weyl semimetals

### Graphene:
- Nanoribbons
- Junctions
- Bilayer systems
- Strain effects
- Edge physics

### Superconductivity:
- Andreev reflection
- Josephson junctions
- Proximity effect
- Majorana fermions
- Topological superconductors

## Best Practices

### System Construction:
- Start simple
- Test with small systems
- Use shape functions
- Validate Hamiltonians
- Visualize systems

### Transport Calculations:
- Energy range planning
- Lead convergence
- Sparse matrix usage
- Parallel energy points
- Error checking

### Performance:
- Sparse matrices
- Appropriate system size
- Lead minimization
- Profiling
- Optimization

## Community and Support
- Open-source (BSD 2-Clause)
- Large community
- Active development
- Mailing list
- Tutorial workshops
- Excellent documentation
- GitHub/GitLab

## Educational Resources
- Comprehensive tutorial
- Example gallery
- Video lectures
- Workshop materials
- API documentation
- Transport theory primer
- Community examples

## Development
- Community project
- Multiple institutions
- Active development
- Regular releases
- Feature additions
- Well-maintained
- Production focus

## Research Impact
Kwant is the standard tool for quantum transport calculations, cited in hundreds of publications across mesoscopic physics, topological materials, and quantum devices.

## Verification & Sources
**Primary sources**:
1. Homepage: https://kwant-project.org/
2. Documentation: https://kwant-project.org/doc/
3. GitLab: https://gitlab.kwant-project.org/kwant/kwant
4. Publications: New J. Phys. 16, 063065 (2014)

**Secondary sources**:
1. Mesoscopic physics literature
2. Quantum transport papers
3. User publications (hundreds)

**Confidence**: VERIFIED - Transport calculation standard

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Repository: ACCESSIBLE
- License: BSD 2-Clause (open-source)
- **Category**: Quantum transport library
- Status: Actively developed
- Community: Large, international
- Specialized strength: Quantum transport calculations, Landauer-Büttiker formalism, flexible geometry construction, scattering matrix, mesoscopic physics, topological transport, superconductivity (BdG), production quality, Python-based, comprehensive documentation, industry standard for transport
