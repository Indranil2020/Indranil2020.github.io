# PyCDT

## Official Resources
- Source Repository: https://github.com/mbkumar/pycdt
- Documentation: https://pycdt.readthedocs.io/
- License: Open source (MIT)

## Overview
**PyCDT** (Python Charged Defect Tools) is a Python package for thermodynamic calculations and error corrections for charged defects in semiconductors and insulators using periodic DFT. It generates inputs for required VASP calculations and processes output to compute defect formation energies with finite-size corrections.

**Scientific domain**: Charged defect formation energy, finite-size corrections, defect thermodynamics  
**Target user community**: Researchers studying point defects in semiconductors and insulators with VASP

## Theoretical Methods
- Defect formation energy calculation
- Potential alignment correction
- Image-charge correction (Makov-Payne, Lany-Zunger)
- Band-filling correction for shallow defects
- Chemical potential determination
- Defect phase diagrams
- Transition level calculation

## Capabilities (CRITICAL)
- Defect formation energy calculation
- Potential alignment correction
- Image-charge correction (Makov-Payne, Lany-Zunger)
- Band-filling correction
- Chemical potential determination
- Defect transition levels
- VASP input generation for defects

**Sources**: GitHub repository, ReadTheDocs, CPC

## Key Strengths

### Comprehensive Corrections:
- Potential alignment
- Image-charge (multiple schemes)
- Band-filling for shallow defects
- All standard finite-size corrections

### VASP Workflow:
- Automatic input generation
- Defect structure creation
- Result processing
- Complete defect workflow

### Thermodynamics:
- Chemical potential phase diagrams
- Defect formation energies
- Transition levels (Epsilon)
- Concentration calculation

## Inputs & Outputs
- **Input formats**:
  - VASP output files
  - Bulk and defect calculations
  - Chemical potential data
  
- **Output data types**:
  - Defect formation energies
  - Transition levels
  - Corrected energies
  - Phase diagrams

## Interfaces & Ecosystem
- **VASP**: Primary DFT backend
- **pymatgen**: Structure handling
- **Python**: Core language

## Performance Characteristics
- **Speed**: Fast (post-processing)
- **Accuracy**: DFT + corrections
- **System size**: Any defect size
- **Memory**: Low

## Computational Cost
- **Analysis**: Seconds
- **VASP pre-requisite**: Hours (separate)
- **Typical**: Efficient

## Limitations & Known Constraints
- **VASP only**: No QE or other code support
- **pymatgen dependency**: Requires pymatgen
- **Supercell approach**: Periodic boundary conditions
- **Limited to standard corrections**: No advanced GW corrections

## Comparison with Other Codes
- **vs doped**: PyCDT is older, doped is newer with more features
- **vs pylada-defects**: PyCDT is VASP-focused, pylada is multi-code
- **vs pymatgen-analysis-defects**: PyCDT is standalone, MP-defects is pymatgen-native
- **Unique strength**: Comprehensive charged defect corrections with multiple image-charge schemes, VASP workflow

## Application Areas

### Semiconductor Defects:
- Point defect formation energies
- Charged defect thermodynamics
- Transition levels
- Defect concentrations

### Photovoltaics:
- Defect tolerance assessment
- Recombination center identification
- Doping efficiency
- Carrier concentration prediction

### Battery Materials:
- Defect-mediated ionic transport
- Vacancy formation energies
- Interstitial stability
- Degradation mechanisms

## Best Practices

### VASP Setup:
- Use sufficiently large supercells
- Include enough k-points for defect cell
- Use same settings for bulk and defect
- Check convergence with supercell size

### Corrections:
- Apply all relevant corrections
- Compare Makov-Payne vs Lany-Zunger
- Check potential alignment convergence
- Validate against known defect levels

## Community and Support
- Open source (MIT)
- ReadTheDocs documentation
- Published in CPC
- Research community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mbkumar/pycdt
2. Documentation: https://pycdt.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- Published methodology: CPC
- Specialized strength: Comprehensive charged defect corrections with multiple image-charge schemes, VASP workflow
