# GASpy (Generalized Adsorption Simulator for Python)

## Official Resources
- Homepage: https://github.com/ulissigroup/GASpy
- Documentation: https://gaspy.readthedocs.io/ (Link broken/Moved)
- Source Repository: https://github.com/ulissigroup/GASpy
- License: MIT License

## Overview
GASpy is a Python package for automating the generation of adsorption structures on catalyst surfaces for high-throughput DFT calculations. Developed by the Ulissi Group (CMU), it focuses on finding the lowest energy adsorption sites for arbitrary adsorbates on arbitrary surfaces using heuristic algorithms and VASP.

**Scientific domain**: Catalysis, high-throughput screening, surface science  
**Target user community**: Catalysis researchers

## Capabilities (CRITICAL)
- **Structure Generation**: Automatically places adsorbates on surfaces.
- **Workflow**: Manages VASP calculations (relaxation).
- **Screening**: Designed for high-throughput screening of catalysts (e.g., CO2 reduction, Nitrogen fixation).
- **Database**: Adds results to a MongoDB database.

**Sources**: GASpy GitHub, Nat. Catal. 1, 663 (2018)

## Inputs & Outputs
- **Input formats**: Bulk structure, adsorbate molecule, miller indices
- **Output data types**: Adsorption energies, structures

## Interfaces & Ecosystem
- **Pymatgen**: Structure handling.
- **FireWorks**: Workflow management.
- **VASP**: DFT engine.

## Workflow and Usage
1. Define bulk and surface (e.g., Cu(111)).
2. Define adsorbate (e.g., CO).
3. GASpy generates initial configurations.
4. Run VASP relaxations.
5. Store lowest energy configuration.

## Performance Characteristics
- Automates the "human intuition" part of finding binding sites.
- Scalable to thousands of surfaces.

## Application Areas
- Catalyst discovery (HEA, intermetallics).
- Adsorption energy databases.

## Community and Support
- Developed by Ulissi Group (CMU/Meta).
- Active research code.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ulissigroup/GASpy
2. Publication: Z. W. Ulissi et al., Nat. Catal. 1, 663 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: MINIMAL
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Adsorption automation
