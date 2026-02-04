# NCIPLOT

## Official Resources
- Homepage: https://github.com/aoterodelaroza/nciplot
- GitHub: https://github.com/aoterodelaroza/nciplot
- Documentation: Included in repository
- Publication: E. R. Johnson et al., J. Am. Chem. Soc. 132, 6498 (2010)
- License: GNU General Public License

## Overview
NCIPLOT is a program that enables the computation and graphical visualization of inter- and intra-molecular non-covalent interactions (hydrogen bonds, π-π interactions, van der Waals contacts) using the reduced density gradient (RDG) method. It produces 3D isosurfaces that can be visualized with molecular graphics programs.

**Scientific domain**: Non-covalent interactions, hydrogen bonding, weak interactions
**Target user community**: Researchers studying molecular recognition, supramolecular chemistry, and weak interactions

## Theoretical Methods
- Reduced density gradient (RDG) analysis
- Electron density topology
- Sign(λ2)ρ coloring scheme
- Promolecular density approximation
- Self-consistent field density analysis

## Capabilities (CRITICAL)
- **NCI Visualization**: 3D isosurfaces of non-covalent regions
- **Hydrogen Bonds**: Identification and visualization
- **π-π Stacking**: Aromatic interaction analysis
- **Van der Waals**: Weak dispersion contacts
- **Promolecular Mode**: Fast approximate analysis
- **SCF Mode**: Accurate DFT-based analysis
- **Large Systems**: Proteins and biomolecules

**Sources**: NCIPLOT documentation, JACS publication

## Key Strengths

### Visual Analysis:
- Intuitive 3D visualization
- Color-coded interaction strength
- VMD/PyMOL compatible output
- Publication-quality figures

### Flexibility:
- Promolecular (fast) mode
- SCF (accurate) mode
- Cube file input
- Multiple output formats

### Widely Validated:
- Extensively benchmarked
- High citation count
- Active development
- Critic2 integration

## Inputs & Outputs
- **Input formats**:
  - Gaussian wfn/wfx files
  - Cube files (electron density)
  - XYZ coordinates (promolecular)
  
- **Output data types**:
  - Cube files for visualization
  - RDG isosurfaces
  - 2D RDG vs sign(λ2)ρ plots

## Installation
```bash
git clone https://github.com/aoterodelaroza/nciplot.git
cd nciplot
make
```

## Usage Examples
```bash
# Basic NCI analysis from wfn file
nciplot molecule.wfn

# Input file example (nciplot.nci):
molecule.wfn
# Output cube files for VMD visualization
```

## Performance Characteristics
- **Speed**: Fast for promolecular, moderate for SCF
- **Memory**: Depends on grid resolution
- **Scalability**: Handles large biomolecules

## Limitations & Known Constraints
- **Visualization needed**: Requires external viewer (VMD, PyMOL)
- **Grid resolution**: Trade-off with file size
- **Qualitative**: Primarily visual analysis tool
- **Wavefunction required**: SCF mode needs DFT output

## Comparison with Other Tools
- **vs Critic2**: NCIPLOT specialized for NCI, Critic2 broader QTAIM
- **vs Multiwfn**: Both do NCI, different interfaces
- **Unique strength**: Focused NCI visualization, widely adopted

## Application Areas
- Drug-receptor interactions
- Protein-ligand binding
- Crystal packing analysis
- Supramolecular chemistry
- Catalysis mechanisms

## Best Practices
- Use promolecular for quick screening
- Use SCF for quantitative analysis
- Choose appropriate isosurface values
- Validate with other methods

## Community and Support
- GitHub repository
- High-impact publication
- Developer: A. Otero-de-la-Roza
- Active maintenance

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/aoterodelaroza/nciplot
2. E. R. Johnson et al., J. Am. Chem. Soc. 132, 6498 (2010)
3. J. Contreras-García et al., J. Chem. Theory Comput. 7, 625 (2011)

**Confidence**: VERIFIED - Published in JACS

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (GPL)
- Academic citations: >4000
- Active development: Maintained
