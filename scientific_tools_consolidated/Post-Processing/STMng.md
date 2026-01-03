# STMng

## Official Resources
- **Homepage**: https://uspex-team.org/en/codes
- **Documentation**: https://uspex-team.org/static/file/Valle2005_Zkrist.pdf
- **Source Repository**: Distributed as part of the USPEX ecosystem or via Mario Valle's website.
- **Developer**: Mario Valle (CSCS)

## Overview
STMng is a visualization tool developed by Mario Valle, specifically designed to interface with **USPEX** (Universal Structure Predictor: Evolutionary Xtallography). It is used for visualizing crystal structures and, notably, the "fingerprints" (Oganov-Valle fingerprints) used in evolutionary algorithms to determine structural similarity. While its name suggests Scanning Tunneling Microscopy (STM) visualization, its primary documented role within the USPEX community is structural analysis and visualization of evolutionary search results.

**Scientific domain**: Crystal structure prediction, Visualization, Structural fingerprints
**Target user community**: USPEX users, crystallographers

## Capabilities
- **Visualization**: Displaying crystal structures from USPEX outputs.
- **Fingerprint Analysis**: Visualizing the Oganov-Valle structural fingerprints (distances in configuration space).
- **Integration**: specialized for the file formats and data structures produced by USPEX.

## Inputs & Outputs
- **Input formats**: USPEX output files (`generations`, specific structure files).
- **Output data types**: 3D visualizations, fingerprint plots.

## Interfaces & Ecosystem
- **USPEX**: The primary parent suite. STMng is listed as a code "well suited for visualizing results of USPEX".

## Verification & Sources
- **Primary Source**: [USPEX Codes Page](https://uspex-team.org/en/codes)
- **Reference**: Valle, M. (2005). *Zeitschrift für Kristallographie*, 220(5-6), 585-588. "Crystal structure visualization and analysis with STMng".
