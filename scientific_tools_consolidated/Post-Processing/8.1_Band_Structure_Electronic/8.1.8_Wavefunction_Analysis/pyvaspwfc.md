# pyvaspwfc

## Official Resources
- Source Repository: https://github.com/liming-liu/pyvaspwfc
- Documentation: Included in repository
- License: Open source

## Overview
**pyvaspwfc** is a Python class for dealing with VASP pseudo-wavefunction file WAVECAR. It can extract planewave coefficients of any Kohn-Sham orbital, perform band unfolding, and visualize wavefunctions in real space via 3D Fourier transform.

**Scientific domain**: VASP WAVECAR analysis, wavefunction visualization, band unfolding  
**Target user community**: Researchers needing to extract and visualize wavefunctions from VASP calculations

## Theoretical Methods
- WAVECAR parsing
- Planewave coefficient extraction
- 3D Fourier transform for real-space wavefunctions
- Band unfolding from supercell wavefunctions
- Pseudo-wavefunction visualization

## Capabilities (CRITICAL)
- WAVECAR parsing and reading
- Planewave coefficient extraction
- Real-space wavefunction visualization
- Band unfolding from supercell
- Charge density calculation
- VASP WAVECAR interface

**Sources**: GitHub repository

## Key Strengths

### WAVECAR Access:
- Direct access to wavefunction data
- Any KS orbital extraction
- Planewave coefficients
- Complete WAVECAR parsing

### Real-Space Visualization:
- 3D Fourier transform
- Real-space wavefunction plots
- Charge density from wavefunctions
- VESTA-compatible output

### Band Unfolding:
- Supercell to primitive cell
- Wavefunction projection
- Spectral weight calculation
- Unfolded band structure

## Inputs & Outputs
- **Input formats**:
  - VASP WAVECAR
  - POSCAR (structure)
  
- **Output data types**:
  - Real-space wavefunctions
  - Planewave coefficients
  - Charge density data
  - Unfolded band data

## Interfaces & Ecosystem
- **VASP**: WAVECAR source
- **NumPy**: Numerical computation
- **Python**: Core language

## Performance Characteristics
- **Speed**: Moderate (WAVECAR is large)
- **Accuracy**: VASP-level
- **System size**: Limited by WAVECAR size
- **Memory**: High (WAVECAR parsing)

## Computational Cost
- **WAVECAR reading**: Seconds to minutes
- **VASP pre-requisite**: Hours (separate)
- **Typical**: Moderate

## Limitations & Known Constraints
- **VASP only**: No QE or other code support
- **WAVECAR required**: Very large file
- **Memory intensive**: Full WAVECAR in memory
- **Gamma-only limitations**: Some WAVECAR types not supported

## Comparison with Other Codes
- **vs pawpyseed**: pyvaspwfc is WAVECAR-focused, pawpyseed is PAW augmentation
- **vs VaspBandUnfolding**: pyvaspwfc has wavefunction viz, VaspBandUnfolding is unfolding only
- **vs VASPBERRY**: pyvaspwfc is wavefunction, VASPBERRY is Berry curvature
- **Unique strength**: WAVECAR parsing with real-space wavefunction visualization and band unfolding

## Application Areas

### Wavefunction Analysis:
- Real-space wavefunction visualization
- Orbital character analysis
- Charge density from wavefunctions
- Bonding analysis

### Band Unfolding:
- Supercell band unfolding
- Spectral weight mapping
- Defect state visualization
- Alloy band structure

### Teaching:
- Wavefunction visualization
- DFT concepts demonstration
- Band structure understanding
- Fourier transform illustration

## Best Practices

### WAVECAR Handling:
- Ensure WAVECAR is complete
- Use appropriate precision
- Check LWAVE flag in VASP
- Manage memory for large systems

### Visualization:
- Use VESTA for 3D visualization
- Choose appropriate isosurface levels
- Compare with charge density
- Validate against known systems

## Community and Support
- Open source on GitHub
- Research code
- Limited documentation
- Example usage provided

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/liming-liu/pyvaspwfc

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Specialized strength: WAVECAR parsing with real-space wavefunction visualization and band unfolding
