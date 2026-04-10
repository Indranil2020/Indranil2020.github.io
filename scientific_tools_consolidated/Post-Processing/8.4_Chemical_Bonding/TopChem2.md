# TopChem2

## Official Resources
- Homepage: https://www.lct.jussieu.fr/pagesperso/pilme/topchempage.html
- Web Interface: http://topchemweb.sorbonne-universite.fr/
- Developer Contact: julien.pilme@sorbonne-universite.fr
- Documentation: Tutorials and manual available from the official page

## Overview
TopChem2 is a standalone quantum chemical topology program for analyzing electron density and related descriptors from wavefunction and cube data. It supports QTAIM analysis together with ELF, NCI, and Fukui-function style descriptors, making it a broad molecular bonding-analysis environment.

**Scientific domain**: QTAIM, ELF, NCI, Fukui descriptors, molecular bonding  
**Target user community**: Quantum chemists, wavefunction analysts, topology users

## Theoretical Methods
- Quantum Theory of Atoms in Molecules (QTAIM)
- Electron Localization Function (ELF)
- Non-Covalent Interaction (NCI) analysis
- Fukui descriptors
- Grid-based numerical topology analysis

## Capabilities (CRITICAL)
- QTAIM analysis from WFN/WFX data
- ELF analysis and basin/topology interpretation
- NCI analysis for intermolecular and intramolecular interactions
- Fukui-descriptor analysis for reactivity interpretation
- Support for wavefunction files and cube-based workflows
- Standalone topology program outside any single electronic-structure code

**Sources**: Official TopChem2 page, TopChemWeb page, tutorial/manual references

## Key Strengths

### Multi-Descriptor Topology:
- QTAIM bonding analysis
- ELF localization analysis
- NCI weak-interaction analysis
- Fukui reactivity descriptors

### Standalone Workflow:
- Works from wavefunction or cube-style data
- Not tied to a single quantum chemistry engine
- Suitable for post-processing studies

### Practical Molecular Analysis:
- Broad descriptor coverage
- Useful for both strong and weak interactions
- Good fit for teaching and research workflows

## Inputs & Outputs
- **Input formats**:
  - WFN and WFX wavefunction files
  - Cube files and related grid data

- **Output data types**:
  - Topological critical point data
  - Basin and descriptor analyses
  - NCI- and ELF-related graphical data
  - Reactivity descriptor outputs

## Workflow and Usage
1. Generate wavefunction or cube files from a quantum chemistry code.
2. Load the data into TopChem2.
3. Select the desired analysis mode: QTAIM, ELF, NCI, or Fukui-related analysis.
4. Inspect topology, interaction regions, and descriptor values.

## Performance Characteristics
- Numerical grid-based post-processing workflow
- Best suited to detailed molecular analysis rather than large-scale automation
- Supports multiple descriptor families in a single environment

## Limitations & Known Constraints
- **Website access**: Official homepage may be intermittently difficult to access from some environments
- **Molecular emphasis**: Primarily a wavefunction and cube post-processing tool
- **Specialized workflow**: Requires prepared WFN/WFX or cube files

## Comparison with Other Tools
- **vs Multiwfn**: Both cover broad post-processing; TopChem2 emphasizes topology-oriented QTAIM/ELF/NCI/Fukui workflows
- **vs Critic2**: TopChem2 is more molecule-focused; Critic2 is stronger for periodic solids
- **Unique strength**: Combines QTAIM, ELF, NCI, and Fukui analysis in one standalone package

## Application Areas
- Molecular bonding analysis
- Weak interaction analysis
- Electron localization studies
- Reactivity interpretation from Fukui descriptors

## Community and Support
- Academic distribution from Sorbonne-associated pages
- Tutorials and manual material available
- Known in quantum-chemical topology workflows

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.lct.jussieu.fr/pagesperso/pilme/topchempage.html
2. Web interface: http://topchemweb.sorbonne-universite.fr/
3. Official page summary describing QTAIM, ELF, NCI, and Fukui analysis

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: KNOWN, but intermittent access from this environment
- Documentation/tutorials: AVAILABLE
- Distribution: AVAILABLE
- Primary use case: Standalone molecular topology analysis across multiple descriptor families
