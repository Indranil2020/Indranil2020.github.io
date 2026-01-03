# pymatgen-db

## Official Resources
- Homepage: https://materialsproject.github.io/pymatgen-db/
- Documentation: https://materialsproject.github.io/pymatgen-db/
- Source Repository: https://github.com/materialsproject/pymatgen-db
- License: Modified BSD License

## Overview
pymatgen-db is a tool for managing MongoDB databases of materials data. It works with Pymatgen to parse calculation results (from VASP, Q-Chem, etc.) and insert them into a database with a structured schema. It allows for powerful querying of materials properties using a Python API or CLI.

**Scientific domain**: Database management, materials informatics  
**Target user community**: Users of FireWorks/Atomate, researchers managing their own data

## Capabilities (CRITICAL)
- **Insertion**: Automates insertion of VASP/Q-Chem outputs into MongoDB.
- **Querying**: Rich query language for filtering by chemistry, properties, and metadata.
- **CLI**: Command-line interface (`mgdb`) for quick database interactions.
- **Schema**: Defines the standard map of JSON keys for calculation data.

**Sources**: pymatgen-db documentation

## Inputs & Outputs
- **Input formats**: Calculation directories, JSON files
- **Output data types**: MongoDB documents

## Interfaces & Ecosystem
- **MongoDB**: Required backend
- **Pymatgen**: Core dependency
- **Atomate**: Uses pymatgen-db logic for database interactions

## Workflow and Usage
1. Configure `db.json`.
2. Insert calculation: `mgdb insert .`
3. Query: `mgdb query --criteria '{"elements": {"$all": ["Li", "O"]}}' --props energy`

## Performance Characteristics
- Efficient management of JSON documents
- Scalable with MongoDB

## Application Areas
- Personal research databases
- Group-level data sharing
- Data backend for high-throughput studies

## Community and Support
- Developed by Materials Project
- Active GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/materialsproject/pymatgen-db

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Materials database management
