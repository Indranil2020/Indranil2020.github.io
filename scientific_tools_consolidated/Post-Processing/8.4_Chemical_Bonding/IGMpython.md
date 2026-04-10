# IGMpython

## Official Resources
- GitHub: https://github.com/bertadenes/IGMpython
- Method basis: Improved implementation of IGMPlot using QM molecular densities
- Related method family: IGM / IGMH analysis

## Overview
IGMpython is a Python implementation of the Independent Gradient Model that uses quantum-mechanical molecular densities provided as cube files. It is designed to reveal and visualize chemical interactions using full and fragment density data, and it produces cube outputs together with a VMD visualization state.

**Scientific domain**: Interaction analysis, IGM post-processing, density-based bonding visualization  
**Target user community**: Computational chemists studying intermolecular and intramolecular interactions from cube data

## Theoretical Methods
- Independent Gradient Model (IGM)
- QM-density-based interaction analysis
- Hessian-eigenvalue coloring inspired by NCI-style visualization
- Fragment-based interaction decomposition from cube densities

## Capabilities (CRITICAL)
- Python 3 implementation using Gaussian cube files as input
- Uses QM molecular densities instead of only promolecular approximations
- Supports fragment-based analysis from full and fragment cube files
- Generates `igm.cub` and coloring cubes such as `mideig.cub` or `diff.cub`
- Automatically generates a VMD state file for visualization

**Sources**: Official GitHub repository and usage documentation

## Key Strengths

### QM-Density Workflow:
- Uses ab initio cube densities
- Fragment-resolved interaction analysis
- Bridges IGM methodology with standard cube workflows

### Practical Outputs:
- Ready-to-visualize cube files
- Automatic VMD state generation
- Multiple coloring options for interpretation

### Lightweight Python Tool:
- Python-based script
- Minimal dependency footprint
- Easy to integrate into analysis workflows

## Inputs & Outputs
- **Input formats**:
  - Full cube files
  - Fragment cube files
  - XYZ files for promolecular-style workflows with `-p`

- **Output data types**:
  - `igm.cub`
  - `mideig.cub` or `diff.cub`
  - VMD visualization state file

## Workflow and Usage
1. Prepare a full cube file and, if needed, fragment cube files on the same grid.
2. Run `IGM.py full.cube -f [fragment.cubes ...]`.
3. Visualize `igm.cub` with the generated coloring cube in VMD.
4. Interpret the interaction regions using the generated surfaces.

## Performance Characteristics
- Lightweight scriptable workflow
- Depends on consistent cube-grid preparation across fragments
- Useful for targeted interaction analysis rather than broad all-in-one post-processing

## Limitations & Known Constraints
- **Grid consistency**: Fragment densities must be represented on the same grid
- **Method scope**: Focused on IGM analysis rather than full topology suites
- **Visualization dependency**: Designed around VMD-oriented output

## Comparison with Other Tools
- **vs IGMPlot**: IGMpython explicitly uses QM molecular densities from cube inputs in a lightweight Python workflow
- **vs NCIPLOT**: Both visualize interaction regions, but IGMpython is centered on the IGM formalism
- **Unique strength**: Simple Python implementation of IGM using full and fragment QM cube densities

## Application Areas
- Weak interaction analysis
- Intramolecular and intermolecular contact visualization
- Fragment-based interaction studies
- VMD-centered interaction mapping

## Community and Support
- Public GitHub repository
- Readme-style installation and usage instructions
- Clearly documented outputs and workflow

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/bertadenes/IGMpython
2. Repository usage documentation describing cube inputs and VMD outputs

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Public repository: ACCESSIBLE
- Installation and usage docs: AVAILABLE
- Primary use case: QM-density-based IGM interaction analysis
