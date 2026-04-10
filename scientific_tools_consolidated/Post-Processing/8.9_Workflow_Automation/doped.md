# doped

## Official Resources
- Source Repository: https://github.com/SMTG-Bham/doped
- Documentation: https://doped.readthedocs.io/
- PyPI: https://pypi.org/project/doped/
- License: MIT License

## Overview
**doped** is a Python software for the generation, simulation, and analysis of defect supercells in solid-state materials. It automates the workflow for point defect calculations using VASP and pymatgen, including symmetry analysis, finite-size corrections, and defect property analysis.

**Scientific domain**: Defect calculations, point defect simulation, defect analysis  
**Target user community**: Researchers studying point defects in semiconductors, insulators, and functional materials

## Theoretical Methods
- Point defect thermodynamics
- Formation energy calculation
- Symmetry analysis of defects
- Finite-size corrections (Kumagai-Oba eFNV, Freysoldt)
- Chemical potential determination
- Fermi level self-consistency
- Charge state analysis
- VASP DFT backend

## Capabilities (CRITICAL)
- Automated defect supercell generation
- Symmetry-inequivalent defect identification
- VASP input file generation
- Competing phase analysis for chemical potentials
- Formation energy calculation
- Finite-size charge correction (eFNV, Freysoldt)
- Defect concentration calculation
- Fermi level vs temperature analysis
- Carrier concentration analysis
- Plotting and analysis tools

**Sources**: GitHub repository, JOSS publication

## Key Strengths

### Automated Workflow:
- End-to-end defect calculation workflow
- Symmetry analysis reduces redundant calculations
- Automatic VASP input generation
- Competing phase determination

### Comprehensive Corrections:
- Freysoldt (FNV) correction
- Kumagai-Oba (eFNV) correction
- Anisotropic dielectric screening
- Band gap alignment

### Rich Analysis:
- Formation energy diagrams
- Concentration vs temperature
- Fermi level determination
- Carrier concentration plots

### Pymatgen Integration:
- Built on pymatgen
- Materials Project API
- Standard file formats
- Extensible framework

## Inputs & Outputs
- **Input formats**:
  - pymatgen Structure objects
  - VASP output files
  - Chemical potential data
  
- **Output data types**:
  - VASP input files for defect calculations
  - Formation energy diagrams
  - Defect concentrations
  - Fermi level vs temperature
  - Charge correction plots

## Interfaces & Ecosystem
- **VASP**: Primary DFT backend
- **pymatgen**: Materials analysis framework
- **Materials Project**: Database access
- **Matplotlib**: Visualization

## Performance Characteristics
- **Speed**: Fast (workflow management)
- **Accuracy**: DFT-level for defect properties
- **System size**: Limited by VASP supercell
- **Automation**: Full workflow automation

## Computational Cost
- **Workflow setup**: Seconds
- **VASP calculations**: Hours to days (separate)
- **Analysis**: Seconds
- **Typical**: Efficient workflow, VASP is bottleneck

## Limitations & Known Constraints
- **VASP only**: No QE or other code support
- **Point defects only**: No extended defects
- **DFT-level**: No GW or hybrid corrections built-in
- **Supercell approach**: Finite-size effects remain

## Comparison with Other Codes
- **vs pymatgen-analysis-defects**: doped is more automated and VASP-focused
- **vs C2DB defect workflow**: doped is general, C2DB is 2D-specific
- **vs PyCDT**: doped is newer with better corrections
- **Unique strength**: End-to-end automated defect calculation workflow with VASP, comprehensive corrections and analysis

## Application Areas

### Semiconductors:
- Native defect formation energies
- Dopant incorporation
- Carrier concentration prediction
- Fermi level engineering

### Photovoltaics:
- Defect tolerance analysis
- Recombination-active defects
- Doping optimization
- Stability assessment

### Battery Materials:
- Cation disorder
- Oxygen vacancy formation
- Electrolyte decomposition
- Degradation mechanisms

### Functional Oxides:
- Oxygen vacancy concentrations
- Redox-active defects
- Conductivity prediction
- Phase stability

## Best Practices

### Supercell Selection:
- Use sufficiently large supercells
- Check convergence with supercell size
- Consider charge state effects
- Validate finite-size corrections

### Chemical Potentials:
- Include all competing phases
- Check stability region
- Use consistent DFT settings
- Consider off-stoichiometry

### Analysis:
- Use appropriate charge corrections
- Consider all charge states
- Include metastable states
- Compare with experiment

## Community and Support
- Open source (MIT License)
- PyPI installation available
- ReadTheDocs documentation
- Developed at University of Birmingham (SMTG)
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/SMTG-Bham/doped
2. Documentation: https://doped.readthedocs.io/
3. M. K. Horton et al., JOSS (2024)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- PyPI: AVAILABLE
- Active development: Ongoing
- Specialized strength: End-to-end automated defect calculation workflow with VASP, comprehensive corrections and analysis
