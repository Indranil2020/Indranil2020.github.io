# FireWorks

## Official Resources
- Homepage: https://materialsproject.github.io/fireworks/
- Documentation: https://materialsproject.github.io/fireworks/
- Source Repository: https://github.com/materialsproject/fireworks
- License: Modified BSD License

## Overview
FireWorks is a workflow software for running, tracking, and managing high-throughput calculations. It is designed to handle complex workflows with dependencies, enabling the automation of large-scale computational tasks. FireWorks uses a centralized database (MongoDB) to store the state of all workflows, allowing for dynamic workflow modification and robust failure recovery.

**Scientific domain**: Workflow management, high-throughput computing
**Target user community**: Computational scientists running large-scale simulations

## Capabilities
- **Centralized Management**: Uses MongoDB to store workflow state, allowing queryable status of all jobs.
- **Dynamic Workflows**: Workflows can modify themselves during execution (e.g., adding new steps based on previous results).
- **Duplicate Detection**: Automatically detects and prevents duplicate runs.
- **Failure Recovery**: Rerun failed jobs, detect "lost" runs.
- **Flexible Execution**: Can run on local machines, supercomputers (SLURM, PBS, etc.), or cloud resources.
- **Python & CLI**: Full control via Python API or command-line interface.

## Interfaces & Ecosystem
- **Part of the Materials Project Stack**: Designed to work with **pymatgen** (analysis), **Custodian** (error correction), and **Atomate** (recipes).
- **Database**: Requires MongoDB.

## Workflow and Usage
A workflow in FireWorks consists of **Fireworks** (individual jobs) connected by dependencies.

1. **Define Workflow**:
   ```python
   from fireworks import Firework, Workflow, LaunchPad
   
   # Define a simple task
   fw1 = Firework(ScriptTask(script="echo 'hello'"), name="hello")
   fw2 = Firework(ScriptTask(script="echo 'world'"), name="world")
   
   # Create workflow: fw2 depends on fw1
   wf = Workflow([fw1, fw2], {fw1: [fw2]}, name="hello_world")
   ```

2. **Add to Database**:
   ```python
   launchpad = LaunchPad.auto_load()
   launchpad.add_wf(wf)
   ```

3. **Execute**:
   Run `rlaunch` (rapid launch) or `qlaunch` (queue launch) on computing resources to pull and execute jobs from the database.

## Application Areas
- High-throughput materials discovery
- Complex simulation pipelines (e.g., relaxation -> static -> band structure)
- Distributed computing tasks

## Verification & Sources
**Primary sources**:
1. Homepage: https://materialsproject.github.io/fireworks/
2. GitHub: https://github.com/materialsproject/fireworks

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Workflow management, Materials Project ecosystem
