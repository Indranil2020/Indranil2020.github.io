# atomate2

## Official Resources
- Homepage: https://materialsproject.github.io/atomate2/
- Documentation: https://materialsproject.github.io/atomate2/
- Source Repository: https://github.com/materialsproject/atomate2
- License: BSD 3-Clause License

## Overview
Atomate2 is the next-generation workflow library for computational materials science, succeeding Atomate. It is built on top of `jobflow` (instead of FireWorks directly) and Pymatgen. It offers a more modern, flexible, and easier-to-use API for defining workflows. It supports VASP, CP2K, ForceField (via LAMMPS), and other codes, with a focus on modularity and dynamic workflow generation.

**Scientific domain**: Workflow automation, high-throughput materials science  
**Target user community**: VASP/CP2K users, next-gen Materials Project contributors

## Capabilities (CRITICAL)
- **Modern Workflow Engine**: Powered by `jobflow`, allowing local execution (without a database) or remote execution (via FireWorks).
- **Dynamic Flows**: Easier creation of dynamic workflows (e.g., loops, conditionals) compared to Atomate 1.
- **Codes Supported**: VASP, CP2K, LAMMPS, Quantum ESPRESSO (via plugins).
- **Output Handling**: Better schema for output documents (JobStore).
- **Interoperability**: Can easily switch between different workflow runners (local, FireWorks, etc.).

**Sources**: Atomate2 documentation

## Inputs & Outputs
- **Input formats**: Pymatgen objects, Jobflow Jobs/Flows
- **Output data types**: JSON/BSON documents (locally or in MongoDB)

## Interfaces & Ecosystem
- **jobflow**: The underlying workflow library
- **pymatgen**: Core materials analysis
- **FireWorks**: Optional backend for remote execution
- **Custodian**: Error handling

## Workflow and Usage
1. Define a job: `job = relax_job(structure)`
2. Run locally (for testing): `run_locally(job)`
3. Or run via FireWorks: `flow = Flow([job]); flow.submit(lpad)`

## Performance Characteristics
- improved serialization compared to Atomate 1
- Faster workflow construction
- Flexible execution models (local vs distributed)

## Application Areas
- High-throughput materials discovery
- Complex simulation pipelines
- Rapid prototyping of workflows

## Community and Support
- Developed by Materials Project team (LBNL)
- Active development on GitHub
- Emerging standard for MP workflows

## Verification & Sources
**Primary sources**:
1. Homepage: https://materialsproject.github.io/atomate2/
2. GitHub: https://github.com/materialsproject/atomate2

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (Materials Project)
- Applications: Next-gen workflows, VASP, CP2K
