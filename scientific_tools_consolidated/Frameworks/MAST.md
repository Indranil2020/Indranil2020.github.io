# MAST (Materials Simulation Toolkit)

## Official Resources
- Homepage: https://github.com/uw-cmg/MAST
- Documentation: https://mast-docs.readthedocs.io/
- Source Repository: https://github.com/uw-cmg/MAST
- License: MIT License

## Overview
MAST (MAterials Simulation Toolkit) is a python tool for managing computational materials science workflows. Developed by the Computational Materials Group at the University of Wisconsin-Madison, it focuses on automating defect workflow calculations, diffusion, and phonon calculations using VASP.

**Scientific domain**: Workflow automation, defects, diffusion  
**Target user community**: VASP users, defect researchers

## Capabilities (CRITICAL)
- **Workflow Management**: Automates dependency handling (e.g., relax -> static -> phonon).
- **Defects**: Specialized workflows for point defect formation energies.
- **Diffusion**: NEB workflows.
- **Recipes**: Input-file driven workflow definition.

**Sources**: MAST documentation, Comp. Mater. Sci. 136, 172 (2017)

## Inputs & Outputs
- **Input formats**: Input recipe files
- **Output data types**: VASP outputs, summary CSVs

## Interfaces & Ecosystem
- **VASP**: Primary engine.
- **Pymatgen**: Used for structure manipulation.

## Workflow and Usage
1. Define a recipe file (text format) specifying ingredients (calculations) and the process (dependencies).
2. Run `mast` command.

## Performance Characteristics
- Lightweight
- Tailored for specific research tasks (defects/diffusion)

## Application Areas
- Defect thermodynamics
- Ion transport in battery materials

## Community and Support
- Developed by Morgan Group (UW-Madison)
- **Status**: Less active than Atomate/AiiDA (Last major updates ~2018)

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/uw-cmg/MAST
2. Publication: T. Mayeshiba et al., Comp. Mater. Sci. 136, 172 (2017)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: MAINTENANCE/INACTIVE
- Applications: Defect workflows, NEB
