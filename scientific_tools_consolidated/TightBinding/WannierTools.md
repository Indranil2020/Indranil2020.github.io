# WannierTools

## Official Resources
- Homepage: http://www.wanniertools.com/
- Documentation: http://www.wanniertools.com/tutorials/
- Source Repository: https://github.com/quanshengwu/wannier_tools
- License: GNU General Public License v3.0

## Overview
WannierTools is a comprehensive software package for investigating topological properties of materials using tight-binding models from Wannier90. Developed by QuanSheng Wu and collaborators, WannierTools calculates topological invariants, surface states, nodal structures, and various topological phenomena. The code has become the standard tool for topological characterization of materials, enabling systematic exploration of topological phases from ab-initio calculations.

**Scientific domain**: Topological materials, Wannier functions, surface states  
**Target user community**: Topological physics researchers, materials scientists, ARPES theorists

## Theoretical Methods
- Topological band theory
- Wilson loop calculations
- Z2 invariants
- Chern numbers
- Mirror Chern numbers
- Weyl/Dirac point detection
- Surface state calculations
- Fermi arc analysis
- Nodal line structures

## Capabilities (CRITICAL)
**Category**: Open-source topological analysis tool
- Topological invariant calculation
- Z2 indices (3D and 2D)
- Chern numbers
- Mirror Chern numbers
- Weyl point finding
- Nodal line detection
- Surface/edge state calculation
- Fermi surface analysis
- Berry curvature
- Anomalous Hall conductivity
- Wannier charge centers
- Landau level spectrum
- ARPES simulation
- Tight-binding from Wannier90
- Production quality

**Sources**: Official website, GitHub, publications

## Key Strengths

### Topological Toolbox:
- Comprehensive invariants
- All major topological classes
- Systematic analysis
- Research and production
- Standard in field

### Wannier90 Integration:
- Direct hr.dat input
- Seamless workflow
- ab-initio to topology
- DFT integration
- Standard pipeline

### Surface States:
- Slab calculations
- Iterative Green's function
- Edge states
- Fermi arcs
- ARPES comparison

### Nodal Structures:
- Weyl/Dirac points
- Nodal lines
- Nodal surfaces
- Systematic search
- Visualization-ready

## Inputs & Outputs
- **Input formats**:
  - wt.in (WannierTools input)
  - wannier90_hr.dat (tight-binding)
  - POSCAR (structure)
  
- **Output data types**:
  - Topological invariants
  - Surface state spectra
  - Band structures
  - Fermi surfaces
  - Berry curvature
  - Gnuplot scripts
  - Visualization data

## Interfaces & Ecosystem

### Wannier90:
- Direct hr.dat input
- Standard workflow
- Tight-binding models
- Universal interface

### Visualization:
- Gnuplot output
- Matplotlib compatible
- VESTA structures
- Publication-ready plots

## Workflow and Usage

### Installation:
```bash
# Clone repository
git clone https://github.com/quanshengwu/wannier_tools.git
cd wannier_tools
# Compile
make
```

### Input File (wt.in):
```fortran
&TB_FILE
Hrfile = 'wannier90_hr.dat'
/

&CONTROL
BulkBand_calc = T
BulkFS_calc = T
BulkGap_cube_calc = T
SlabBand_calc = T
Z2_3D_calc = T
WeylPoints_calc = T
/

&SYSTEM
NumOccupied = 18
SOC = 1
E_FERMI = 0.0
/

&PARAMETERS
Nk1 = 101
Nk2 = 101
Nk3 = 101
NP = 2
Gap_threshold = 0.01
/

KPATH_BULK
4
G 0.0 0.0 0.0 X 0.5 0.0 0.0
X 0.5 0.0 0.0 M 0.5 0.5 0.0
M 0.5 0.5 0.0 G 0.0 0.0 0.0
G 0.0 0.0 0.0 Z 0.0 0.0 0.5

SURFACE
1 0 0
/

KPATH_SLAB
2
K 0.33 0.67 G 0.0 0.0
G 0.0 0.0 M 0.5 0.5
```

### Run WannierTools:
```bash
# Execute
wt.x
```

### Visualize Results:
```bash
# Plot bulk bands
gnuplot bulkek.gnu

# Plot surface states
gnuplot surfdos_l.gnu

# Plot Fermi surface
gnuplot fs.gnu
```

## Advanced Features

### Z2 Invariants:
- 3D topological insulators
- 2D topological insulators
- Wilson loop method
- Four Z2 indices (ν0;ν1ν2ν3)
- Automated calculation

### Weyl Physics:
- Weyl point detection
- Chirality calculation
- Fermi arc connection
- Type-I and Type-II Weyl
- Systematic search

### Berry Curvature:
- Momentum-space distribution
- Integration for Chern number
- Anomalous Hall conductivity
- Berry dipole
- Visualization

### Landau Levels:
- Magnetic field response
- Landau fan diagram
- Quantum oscillations
- Topological signatures

## Performance Characteristics
- **Speed**: Fast (post-Wannier90)
- **Accuracy**: Tight-binding quality
- **System size**: Any (uses TB model)
- **Purpose**: Topological analysis
- **Typical**: Minutes to hours

## Computational Cost
- Post-Wannier90 processing
- k-point mesh dependent
- Surface states most expensive
- Efficient algorithms
- Production capable

## Limitations & Known Constraints
- **Requires Wannier90**: Not standalone DFT
- **Tight-binding**: Quality depends on MLWFs
- **Surface states**: Computational cost for large slabs
- **k-mesh**: Dense grids needed
- **Interpretation**: Requires physics knowledge

## Comparison with Other Tools
- **vs Z2Pack**: WannierTools comprehensive, Z2Pack specialized
- **vs WannierBerri**: WannierTools topology, WannierBerri transport
- **Unique strength**: Most comprehensive topological toolbox, standard for topological analysis, Wannier90 integration

## Application Areas

### Topological Insulators:
- Z2 classification
- Surface states
- Edge states
- 2D and 3D TIs
- Material prediction

### Topological Semimetals:
- Weyl semimetals
- Dirac semimetals
- Nodal line semimetals
- Type-II Weyl
- Fermi arcs

### Topological Characterization:
- Material screening
- Phase classification
- Symmetry analysis
- Database generation
- High-throughput

### ARPES Theory:
- Surface state prediction
- Comparison with experiment
- Spectral functions
- Momentum-space features

## Best Practices

### Wannier Functions:
- Quality MLWFs from Wannier90
- Appropriate projections
- Converged tight-binding
- Validated band structure

### k-Point Grids:
- Dense for topology
- Convergence testing
- Surface states need more
- Balance accuracy/cost

### Topological Analysis:
- Check multiple invariants
- Symmetry considerations
- Gap requirements
- Physical interpretation

## Community and Support
- Open-source (GPL v3)
- Active development
- GitHub repository
- User manual
- Example gallery
- Publications
- Growing community

## Educational Resources
- Comprehensive tutorials
- Example inputs
- Gallery of topological materials
- Publication list
- Topological theory background
- Visualization examples

## Development
- QuanSheng Wu (lead, ETH Zurich)
- Alexey Soluyanov group
- Active development
- Regular updates
- Feature additions
- Community contributions

## Research Impact
WannierTools has become the standard tool for topological analysis of materials from first principles, cited in hundreds of publications on topological insulators, Weyl semimetals, and other topological phases.

## Verification & Sources
**Primary sources**:
1. Homepage: http://www.wanniertools.com/
2. GitHub: https://github.com/quanshengwu/wannier_tools
3. Publications: Comp. Phys. Comm. 224, 405 (2018)

**Secondary sources**:
1. Topological materials papers
2. User publications
3. ARPES literature

**Confidence**: CONFIRMED - Standard topological tool

**Verification status**: ✅ CONFIRMED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: GPL v3 (open-source)
- **Category**: Open-source topological analysis tool
- Status: Actively developed
- Institution: ETH Zurich (Soluyanov group)
- Specialized strength: Comprehensive topological invariant calculations, Z2/Chern numbers, Weyl point detection, surface states, Fermi arcs, Wannier90 integration, standard for topological materials analysis, ARPES simulation, nodal structures, production quality, visualization-ready output
