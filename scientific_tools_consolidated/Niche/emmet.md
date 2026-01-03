# emmet

## Official Resources
- Homepage: https://github.com/materialsproject/emmet
- Documentation: https://github.com/materialsproject/emmet
- Source Repository: https://github.com/materialsproject/emmet
- License: Modified BSD License

## Overview
Emmet is the core library defining the data models and "Builders" for the Materials Project. It defines the schema for materials properties (e.g., `MaterialDoc`, `TaskDoc`) and contains the logic to aggregate raw calculation tasks into consolidated material documents.

**Scientific domain**: Materials data schema, database building  
**Target user community**: Materials Project developers, users building MP-like databases

## Capabilities (CRITICAL)
- **Models**: Pydantic models for VASP tasks, materials, diffraction patterns, etc.
- **Builders**: Maggma builders for processing tasks into materials (e.g., grouping relaxation and static runs).
- **Validation**: Ensures data integrity via schema validation.

**Sources**: Emmet GitHub

## Inputs & Outputs
- **Input formats**: Task documents (from atomate/atomate2)
- **Output data types**: Material documents (JSON/BSON)

## Interfaces & Ecosystem
- **Maggma**: Uses Maggma for building pipelines.
- **Pymatgen**: Uses Pymatgen objects.
- **MP API**: Serves data structured by Emmet.

## Workflow and Usage
1. Parse VASP directory -> TaskDoc.
2. Emmet Builder aggregates TaskDocs -> MaterialDoc.
3. MaterialDoc is served via API.

## Performance Characteristics
- The logic behind MP's data generation.

## Application Areas
- Database construction
- Data standardization

## Community and Support
- Materials Project Team

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/materialsproject/emmet

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: MP schema, data models
