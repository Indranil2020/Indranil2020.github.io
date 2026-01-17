# SAX (Likely SAXS - Experimental Technique)

## Official Resources
- Homepage: N/A (Experimental Technique)
- Documentation: Standard Crystallography/Spectroscopy Textbooks
- Source Repository: N/A
- License: N/A

## Overview
Based on extensive verification, "SAX" does not appear to be a standalone Time-Dependent Density Functional Theory (TDDFT) or excited-state software code. It is highly probable that this entry is a misclassification of the experimental technique **Small-Angle X-ray Scattering (SAXS)**, or a typo for the **SAXIS** parameter in VASP (Spin-Quantization Axis).

**Scientific domain**: Experimental Characterization (SAXS), not Computational Chemistry Software
**Status**: **METHOD / MISIDENTIFIED**

## Analysis of Potential Meanings

### 1. Small-Angle X-ray Scattering (SAXS)
- **Context**: A widely used experimental technique for structural characterization of non-crystalline systems (proteins, polymers, colloids) at the nanoscale.
- **Relevance**: Often cited alongside computational spectroscopy (to compare calculated vs experimental scattering profiles), which may have led to its inclusion in a list of "tools".
- **Software**: There are software codes *for* SAXS analysis (e.g., ATSAS, FoXS, CRYSOL), but "SAX" itself is the method.

### 2. SAXIS (VASP Parameter)
- **Context**: `SAXIS` is a specific input tag in the VASP code to define the direction of the global spin-quantization axis for non-collinear magnetic calculations.
- **Relevance**: Appears in input files for advanced DFT calculations.

### 3. Typo for SALMON
- **Context**: SALMON is a real-space TDDFT code. "SAX" shares phonetic similarity or could be a truncation.

## Recommendation
This entry should likely be removed from a list of *software codes* or re-labeled as a "Method" if the list intends to cover techniques. If you are looking for software to simulate SAXS patterns from atomic structures:
- **CRYSOL/ATSAS**: For biological macromolecules
- **FoXS**: Fast X-ray Scattering profile computation
- **CP2K/Quantum ESPRESSO**: Can calculate structure factors which relate to scattering

## Verification & Sources
**Primary sources**:
1. Literature search for "SAX code TDDFT" -> Yields zero software results.
2. Literature search for "SAX technique" -> Yields Small-Angle X-ray Scattering.
3. VASP Wiki (SAXIS tag).

**Confidence**: VERIFIED (as non-software)

**Verification status**: ❌ NOT A SOFTWARE
- **Correction**: Entry corresponds to an experimental method or input parameter, not a standalone simulation package.
