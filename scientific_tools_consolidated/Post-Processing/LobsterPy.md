# LobsterPy

## Official Resources
- Homepage: https://pypi.org/project/lobsterpy/
- Documentation: https://lobsterpy.readthedocs.io/en/latest/
- Source Repository: https://github.com/JaGeo/LobsterPy
- License: MIT License

## Overview
LobsterPy is a Python package for analyzing and plotting data from the LOBSTER code (Local Orbital Basis Suite Towards Electronic-Structure Reconstruction). It automates the analysis of chemical bonding information (COHP, COOP) and provides tools for visualizing bonding properties, stability analysis, and extracting key bonding descriptors automatically.

**Scientific domain**: Chemical bonding analysis, automated plotting, electronic structure  
**Target user community**: Chemists, materials scientists using LOBSTER

## Capabilities (CRITICAL)
- **Automated Plotting**: Generation of publication-quality COHP/COOP/DOS plots
- **Analysis**: Automatic identification of relevant bonds and interactions
- **Descriptions**: Automatic generation of text descriptions of bonding situations
- **Integration**: Works with pymatgen objects
- **Batch Processing**: Tools for high-throughput bonding analysis
- **Machine Learning**: Features for generating bonding descriptors for ML

**Sources**: LobsterPy documentation, GitHub repository

## Inputs & Outputs
- **Input formats**: LOBSTER output files (COHPCAR.lobster, ICOHPLIST.lobster, DOSCAR.lobster)
- **Output data types**: Matplotlib figures, JSON summaries, Python objects

## Interfaces & Ecosystem
- **LOBSTER**: Primary backend code
- **Pymatgen**: Integration for structure and electronic data
- **Python**: Library and command-line tool (`lobsterpy`)

## Workflow and Usage
1. Run LOBSTER calculation.
2. Run `lobsterpy description` to get an automated text summary of bonding.
3. Run `lobsterpy plot` to generate COHP plots for relevant bonds.
4. Use Python API for custom analysis.

## Performance Characteristics
- Fast Python-based post-processing
- Efficient handling of large LOBSTER output files

## Application Areas
- High-throughput bonding analysis
- Understanding crystal stability
- Educational visualization of bonding

## Community and Support
- Open-source (MIT)
- Developed by Janine George group (BAM/Jena)
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/JaGeo/LobsterPy
2. Documentation: https://lobsterpy.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (George Group)
- Applications: Automated COHP analysis, LOBSTER post-processing

