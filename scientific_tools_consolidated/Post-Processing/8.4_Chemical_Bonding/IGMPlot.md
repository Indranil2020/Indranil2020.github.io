# IGMPlot

## Official Resources
- Homepage: http://igmplot.univ-reims.fr/
- Publication: J. Hénon et al., J. Comput. Chem. 44, 1938 (2023)
- Method background: IGM and IGMH literature referenced from the official site

## Overview
IGMPlot is a program for identifying, characterizing, and quantifying molecular interactions using the Independent Gradient Model. It supports analysis of interactions ranging from weak non-covalent contacts to stronger covalent or coordinative interactions, with graphical interpretation based on promolecular or quantum-mechanical density.

**Scientific domain**: Interaction analysis, non-covalent interactions, bonding visualization  
**Target user community**: Computational chemists studying intermolecular and intramolecular interactions

## Theoretical Methods
- Independent Gradient Model (IGM)
- Independent Gradient Model based on Hirshfeld partition (IGMH)
- Interaction region visualization via δg-type descriptors
- Critical-point assisted interaction interpretation

## Capabilities (CRITICAL)
- Detection and visualization of weak and strong interactions
- Analysis of intramolecular and intermolecular contacts
- Support for promolecular-density and QM-density workflows
- Quantification of interaction signatures over a broad bonding range
- Useful analysis of hydrogen bonds, dispersion, steric effects, and coordination interactions
- Modern published implementation dedicated to interaction analysis

**Sources**: Official IGMPlot site and JCC publication

## Key Strengths

### Broad Interaction Coverage:
- Weak non-covalent contacts
- Hydrogen bonding
- Dispersion-driven contacts
- Metal coordination and stronger bonding regions

### Practical Visualization:
- Interaction maps
- Quantitative descriptors
- Supports interpretation across many chemical systems
- Designed specifically for interaction analysis workflows

### Modern Methodology:
- Published recent implementation
- Builds on IGM and IGMH developments
- Useful complement to NCI-focused tools

## Inputs & Outputs
- **Input formats**:
  - Promolecular or QM-density derived inputs as supported by the official program workflow

- **Output data types**:
  - Interaction-region maps
  - Quantification metrics for interactions
  - Graphical files for visualization and analysis

## Workflow and Usage
1. Prepare the required molecular density-related input.
2. Run IGMPlot for IGM or IGMH analysis.
3. Generate interaction maps and descriptors.
4. Interpret weak and strong interactions using the plotted regions and metrics.

## Performance Characteristics
- Specialized post-processing for interaction analysis
- Intended for interpretation and visualization rather than full general-purpose wavefunction processing
- Modern tool focused on user-oriented interaction mapping

## Limitations & Known Constraints
- **Method scope**: Focused on IGM/IGMH interaction analysis rather than full QTAIM basin analysis
- **Workflow specialization**: Best used when interaction visualization is the main goal
- **Complementary role**: Often used alongside QTAIM, NCI, or orbital analyses

## Comparison with Other Tools
- **vs NCIPLOT**: Both analyze interactions; IGMPlot emphasizes IGM/IGMH descriptors and broad interaction quantification
- **vs TopChem2**: IGMPlot is more specialized for interaction mapping, while TopChem2 spans multiple topology descriptors
- **Unique strength**: Dedicated IGM/IGMH program for weak-to-strong interaction analysis

## Application Areas
- Non-covalent interaction analysis
- Intramolecular contact analysis
- Metal coordination studies
- Bonding visualization in complex molecular systems

## Community and Support
- Official project website
- Recent peer-reviewed JCC publication
- Research-oriented software for interaction analysis

## Verification & Sources
**Primary sources**:
1. Homepage: http://igmplot.univ-reims.fr/
2. J. Hénon et al., J. Comput. Chem. 44, 1938 (2023)
3. Official program description for IGM/IGMH analysis

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Publication: AVAILABLE
- Distribution: AVAILABLE
- Primary use case: IGM/IGMH-based interaction analysis and visualization
