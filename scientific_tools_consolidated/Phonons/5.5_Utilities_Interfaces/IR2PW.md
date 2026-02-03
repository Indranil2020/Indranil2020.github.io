# IR2PW

## Official Resources
- Homepage: https://github.com/zjwang11/IR2PW
- Source Repository: https://github.com/zjwang11/IR2PW
- License: Open Source

## Overview
IR2PW (Irreducible Representations to Plane Wave) is a tool for analyzing irreducible representations of phonon modes from plane-wave DFT calculations. It provides symmetry analysis and mode classification for phonon calculations.

**Scientific domain**: Phonon symmetry, irreducible representations  
**Target user community**: Researchers analyzing phonon symmetry from DFT

## Theoretical Methods
- Irreducible representation analysis
- Plane-wave phonon symmetry
- Group theory
- Mode classification
- Character analysis

## Capabilities (CRITICAL)
- Phonon irrep analysis
- Plane-wave code interface
- Symmetry classification
- Mode labeling
- VASP compatibility

## Key Strengths

### Symmetry Analysis:
- Irrep determination
- Mode classification
- Selection rules
- Group theory based

### DFT Integration:
- VASP interface
- Plane-wave codes
- Standard workflow

## Inputs & Outputs
- **Input formats**:
  - VASP phonon output
  - OUTCAR files
  - Structure files
  
- **Output data types**:
  - Irrep labels
  - Mode symmetries
  - Classification tables

## Interfaces & Ecosystem
- **VASP**: Primary interface
- **Python**: Implementation


## Advanced Features
- **Irreducible representation analysis**: Complete symmetry classification
- **Character table analysis**: Group theory-based mode assignment
- **Selection rules**: Raman and IR activity determination
- **VASP integration**: Direct parsing of VASP output
- **Mode labeling**: Automatic symmetry labels

## Performance Characteristics
- Post-processing tool: Fast
- Depends on phonon calculation size
- Python-based implementation

## Computational Cost
- Phonon calculation: External (VASP)
- IR2PW analysis: Fast (seconds to minutes)
- Overall: Minimal overhead

## Best Practices
- Ensure proper symmetry in DFT calculation
- Use appropriate k-point sampling
- Validate against known symmetry assignments
- Compare with experimental spectroscopy

## Limitations & Known Constraints
- VASP-specific
- Limited documentation
- Expert tool
- Specific use case

## Application Areas
- Phonon symmetry analysis
- Mode assignment
- Selection rules
- Spectroscopy interpretation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/zjwang11/IR2PW

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
