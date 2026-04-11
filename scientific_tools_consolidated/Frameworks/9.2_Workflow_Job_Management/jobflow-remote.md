# jobflow-remote

## Official Resources
- Homepage: https://github.com/Matgenix/jobflow-remote
- Documentation: https://matgenix.github.io/jobflow-remote/
- Source Repository: https://github.com/Matgenix/jobflow-remote
- License: Modified BSD License

## Overview
Jobflow-remote (formerly `jobflow-runners` or similar concepts) is a runner for `jobflow` that enables execution on remote resources (like SLURM clusters) without the full complexity of FireWorks' database-polling model. It is a lightweight alternative for submitting jobflow Flows to HPC queues.

**Scientific domain**: Remote execution, HPC job submission  
**Target user community**: Jobflow users needing HPC execution

## Capabilities (CRITICAL)
- **Remote Submission**: Submit jobflow Jobs to SLURM/PBS queues.
- **Lightweight**: Does not require a central MongoDB server for job management (unlike FireWorks).
- **Job Management**: Tracks status of submitted jobs.
- **Integration**: Works seamlessly with `atomate2` and `jobflow`.

**Sources**: jobflow-remote documentation

## Inputs & Outputs
- **Input formats**: Jobflow Flows
- **Output data types**: Job results (files/database)

## Interfaces & Ecosystem
- **jobflow**: The parent workflow library
- **HPC Schedulers**: SLURM, PBS, etc.

## Workflow and Usage
1. Configure remote host (SSH, scheduler type).
2. `submit_flow(flow, worker=remote_worker)`

## Performance Characteristics
- Low latency submission
- Simplified architecture for smaller groups or individual users

## Application Areas
- Running atomate2 workflows on clusters
- Managing calculations from a local machine

## Community and Support
- Developed by Matgenix (Hautier group) and contributors
- Open source

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/Matgenix/jobflow-remote

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN
- Development: ACTIVE
- Applications: Remote execution for jobflow
