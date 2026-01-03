# Wannier90

## Official Resources
- Homepage: https://wannier.org/
- Documentation: https://wannier.org/documentation/
- Source Repository: https://github.com/wannier-developers/wannier90
- License: GNU General Public License v2.0

## Overview
Wannier90 is the community standard code for computing maximally-localized Wannier functions (MLWFs) from electronic structure calculations. Developed by an international collaboration, Wannier90 transforms delocalized Bloch states from DFT calculations into localized Wannier representations, enabling tight-binding model construction, topological analysis, transport calculations, and many other applications. The code interfaces with virtually all major DFT packages and is essential for modern electronic structure workflows.

**Scientific domain**: Wannier functions, tight-binding models, electronic structure  
**Target user community**: DFT users, materials scientists, topological physics, transport calculations

## Theoretical Methods
- Maximally-localized Wannier functions (MLWFs)
- Disentanglement algorithms
- Wannier interpolation
- Berry phase calculations
- Topological invariants
- Tight-binding Hamiltonian construction
- ab-initio tight-binding

## Capabilities (CRITICAL)
**Category**: Open-source Wannier function code
- MLWF construction
- Disentanglement procedures
- Wannier interpolation
- Band structure interpolation
- Density of states
- Berry curvature
- Anomalous Hall conductivity
- Orbital magnetization
- Topological invariants (Z2, Chern)
- Interface with major DFT codes
- Parallelization (MPI)
- Post-processing utilities
- Production quality
- Community standard

**Sources**: Official website, documentation, publications

## Key Strengths

### Community Standard:
- Most widely used
- DFT integration essential
- Thousands of citations
- Production quality
- Extensive validation

### Universal DFT Interface:
- VASP, Quantum ESPRESSO
- ABINIT, Siesta, CASTEP
- Wien2k, FLEUR, Elk
- CP2K, OpenMX
- Universal tool

### Comprehensive Features:
- Band interpolation
- Topological properties
- Transport calculations
- Spectroscopy
- Materials analysis

### Active Development:
- International collaboration
- Regular releases
- New features
- Community support
- Best practices

## Inputs & Outputs
- **Input formats**:
  - .win (Wannier90 input)
  - DFT overlap matrices (from DFT codes)
  - MLWF projections
  - k-point grids
  
- **Output data types**:
  - Wannier functions
  - Tight-binding Hamiltonians (.tb files)
  - Interpolated band structures
  - Berry curvature
  - hr.dat (real-space Hamiltonian)
  - Topological invariants

## Interfaces & Ecosystem

### DFT Codes:
- VASP (vasp2wannier90)
- Quantum ESPRESSO (pw2wannier90)
- ABINIT (native)
- CASTEP (native)
- Siesta, Wien2k, FLEUR
- Nearly universal

### Post-Processing:
- WannierTools (topological)
- WannierBerri (Berry properties)
- BoltzWann (transport)
- BerryPI (polarization)
- Many downstream tools

## Workflow and Usage

### Installation:
```bash
# Download from wannier.org or GitHub
git clone https://github.com/wannier-developers/wannier90.git
cd wannier90
make
```

### DFT Calculation (Quantum ESPRESSO):
```bash
# SCF calculation
pw.x < scf.in > scf.out

# NSCF on fine k-grid
pw.x < nscf.in > nscf.out

# Generate inputs for Wannier90
wannier90.x -pp silicon

# Compute overlaps
pw2wannier90.x < pw2wan.in > pw2wan.out
```

### Wannier90 Input (silicon.win):
```
num_wann = 4
num_iter = 100

begin projections
Si:sp3
end projections

begin unit_cell_cart
bohr
0.0  5.13  5.13
5.13 0.0   5.13
5.13 5.13  0.0
end unit_cell_cart

begin atoms_frac
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25
end atoms_frac

mp_grid = 4 4 4

begin kpoints
# k-point list
end kpoints

bands_plot = true
write_hr = true
```

### Run Wannier90:
```bash
# Wannierize
wannier90.x silicon
```

### Interpolate Band Structure:
```bash
# Add to .win file:
bands_plot = true
begin kpoint_path
L 0.5 0.5 0.5 G 0.0 0.0 0.0
G 0.0 0.0 0.0 X 0.5 0.0 0.5
end kpoint_path

# Rerun
wannier90.x silicon

# Plot bands
gnuplot silicon_band.gnu
```

## Advanced Features

### Disentanglement:
- Energy window selection
- Inner/outer windows
- Projection optimization
- Subspace selection
- Quality control

### Berry Phase Properties:
- Berry curvature
- Anomalous Hall conductivity
- Orbital magnetization
- Chern numbers
- Modern theory of polarization

### Topological Invariants:
- Z2 invariants
- Chern numbers
- Mirror Chern numbers
- Weyl points
- Topological characterization

### Transport:
- Boltzmann transport (with BoltzWann)
- Conductivity tensors
- Seebeck coefficient
- Electronic structure for transport

## Performance Characteristics
- **Speed**: Efficient, post-DFT
- **Accuracy**: High-quality MLWFs
- **System size**: Any (post-processing)
- **Purpose**: Production standard
- **Typical**: Minutes to hours post-DFT

## Computational Cost
- Post-processing (after DFT)
- DFT calculation dominant
- Wannierization fast
- k-point dependent
- Production efficient

## Limitations & Known Constraints
- **Requires DFT**: Not standalone
- **Disentanglement**: Can be tricky
- **Projection choice**: User expertise
- **Localization**: Not always achievable
- **Gauge freedom**: Choices matter

## Comparison with Other Tools
- **Community standard**: No real alternative for MLWFs
- **Complements**: DFT packages
- **Downstream**: WannierTools, WannierBerri, etc. build on it
- **Unique position**: Essential post-DFT tool

## Application Areas

### Materials Science:
- Electronic structure
- Band structures
- Tight-binding models
- Properties calculations
- Universal application

### Topological Materials:
- Topological insulators
- Weyl semimetals
- Chern insulators
- Z2 invariants
- Berry curvature

### Transport:
- Boltzmann transport
- Thermoelectrics
- Hall effects
- Conductivity
- Device modeling

### Spectroscopy:
- Optical properties
- ARPES simulation
- Photoemission
- Interband transitions

## Best Practices

### Projection Choice:
- Physical atomic orbitals
- Symmetry-adapted
- Test different projections
- Validate localization

### Disentanglement:
- Appropriate energy windows
- Sufficient k-points
- Check spreads
- Convergence testing

### k-Point Grid:
- Dense for interpolation
- Converge with grid
- Symmetry exploitation
- Computational balance

## Community and Support
- Open-source (GPL v2)
- International collaboration
- Large user community
- Active mailing list
- Workshops and schools
- Extensive documentation
- GitHub repository

## Educational Resources
- Comprehensive tutorial
- User guide
- Example inputs
- Workshop materials
- Scientific papers
- Video lectures
- MLWF theory primers

## Development
- International collaboration
- Multiple institutions
- Active development
- Regular releases (v3.x)
- Feature additions
- Community-driven
- Best-practice standard

## Research Impact
Wannier90 is cited in thousands of publications and is essential infrastructure for modern electronic structure calculations, enabling topological analysis, transport properties, and tight-binding model construction across materials science and condensed matter physics.

## Verification & Sources
**Primary sources**:
1. Homepage: https://wannier.org/
2. Documentation: https://wannier.org/documentation/
3. GitHub: https://github.com/wannier-developers/wannier90
4. Publications: Comp. Phys. Comm. 178, 685 (2008); 185, 2309 (2014)

**Secondary sources**:
1. User publications (thousands)
2. DFT package documentation
3. Topological materials literature
4. Transport calculations

**Confidence**: CONFIRMED - Community standard

**Verification status**: ✅ CONFIRMED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: GPL v2 (open-source)
- **Category**: Open-source Wannier function code
- Status: Actively developed
- Community: Very large, international
- Specialized strength: Community standard for maximally-localized Wannier functions, universal DFT interface, tight-binding model construction, topological invariants, Berry phase properties, essential post-DFT tool, production quality, comprehensive features, thousands of citations, enables downstream tools ecosystem
