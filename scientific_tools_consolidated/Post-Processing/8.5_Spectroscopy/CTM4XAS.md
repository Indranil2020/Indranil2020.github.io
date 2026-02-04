# CTM4XAS

## Official Resources
- Homepage: https://anorg.chem.uu.nl/CTM4XAS/
- Download: https://anorg.chem.uu.nl/CTM4XAS/software.html
- Publication: E. Stavitski, F.M.F. de Groot, Micron 41, 687 (2010)
- License: Academic use

## Overview
CTM4XAS (Charge Transfer Multiplet for X-ray Absorption Spectroscopy) is a program for simulating L-edge and M-edge X-ray absorption spectra of transition metal compounds using atomic multiplet theory and charge transfer effects. It provides a user-friendly interface for calculating XAS and EELS spectra.

**Scientific domain**: XAS multiplet simulation, L-edge spectroscopy
**Target user community**: X-ray spectroscopists studying transition metals

## Theoretical Methods
- Atomic multiplet theory
- Crystal field effects
- Charge transfer multiplets
- Spin-orbit coupling
- 2p-3d and 3d-4d transitions
- Ligand field theory

## Capabilities (CRITICAL)
- **L-edge XAS**: 2p → 3d transitions
- **M-edge XAS**: 3p → 3d, 3d → 4f transitions
- **EELS**: Electron energy loss spectra
- **Crystal Field**: Various symmetries
- **Charge Transfer**: LMCT, MLCT effects
- **GUI Interface**: User-friendly operation
- **Fitting**: Spectral fitting capabilities

**Sources**: CTM4XAS documentation, Micron publication

## Key Strengths

### Multiplet Theory:
- Full multiplet treatment
- Charge transfer effects
- Crystal field splitting
- Validated methodology

### User-Friendly:
- GUI interface
- Parameter input
- Quick calculations
- Visual output

### Transition Metals:
- 3d transition metals
- 4d and 5d elements
- Rare earth edges
- Various oxidation states

## Inputs & Outputs
- **Input formats**:
  - GUI parameter entry
  - Crystal field parameters
  - Slater integrals
  
- **Output data types**:
  - XAS spectra
  - EELS spectra
  - Energy level diagrams
  - ASCII data files

## Installation
```
# Download from CTM4XAS website
# Windows executable available
# Academic license required
```

## Usage Examples
```
# Typical workflow:
1. Select element and edge (e.g., Fe L2,3)
2. Choose symmetry (e.g., Oh, D4h)
3. Set crystal field parameters (10Dq)
4. Adjust charge transfer parameters if needed
5. Calculate spectrum
6. Compare with experiment
```

## Performance Characteristics
- **Speed**: Fast multiplet calculations
- **Memory**: Minimal requirements
- **Ease of use**: GUI-based

## Limitations & Known Constraints
- **Semi-empirical**: Requires parameter fitting
- **Windows only**: Limited platform support
- **Atomic approach**: Solid-state effects approximate
- **Documentation**: Could be more extensive

## Comparison with Other Tools
- **vs EDRIXS**: CTM4XAS GUI, EDRIXS Python/exact
- **vs Quanty/Crispy**: Different implementations
- **vs FEFF**: CTM4XAS multiplet, FEFF real-space MS
- **Unique strength**: User-friendly L-edge XAS

## Application Areas
- Transition metal oxides
- Coordination chemistry
- Catalysis research
- Magnetic materials
- Battery materials

## Best Practices
- Start with literature parameters
- Validate against known compounds
- Consider charge transfer effects
- Compare multiple symmetries

## Community and Support
- Utrecht University development
- Academic distribution
- Published methodology
- F.M.F. de Groot (developer)

## Verification & Sources
**Primary sources**:
1. Homepage: https://anorg.chem.uu.nl/CTM4XAS/
2. E. Stavitski, F.M.F. de Groot, Micron 41, 687 (2010)

**Confidence**: VERIFIED - Standard XAS tool

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Publication: Micron (2010)
- Academic citations: >1000
- Method: Charge transfer multiplet
