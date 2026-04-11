# MyQueue

## Official Resources
- Homepage: https://myqueue.readthedocs.io/
- Documentation: https://myqueue.readthedocs.io/
- Source Repository: https://gitlab.com/schmidt-group/myqueue
- License: GNU General Public License v3.0

## Overview
MyQueue is a lightweight, frontend-agnostic task scheduler and workflow manager. It is designed to manage high-throughput calculations on local machines or HPC clusters (SLURM, PBS, etc.) with a simple command-line interface. Unlike complex workflow engines, MyQueue focuses on simplicity, using a folder-based structure where folders represent tasks.

**Scientific domain**: Job scheduling, simple workflows, high-throughput  
**Target user community**: ASE/GPAW users, computational physicists

## Capabilities (CRITICAL)
- **Folder-based Tasks**: A folder with a script (e.g., `calculate.py`) is a task.
- **Dependency Tracking**: Dependencies defined by folder hierarchy or explicit links.
- **Scheduler Agnostic**: Unified interface for SLURM, PBS, LSF, and local execution.
- **Command Line**: Simple commands (`mq submit`, `mq list`, `mq kick`).
- **Python Integration**: Works well with ASE scripts.
- **Restart**: Easy restart of failed tasks.

**Sources**: MyQueue documentation

## Inputs & Outputs
- **Input formats**: Python scripts inside directories
- **Output data types**: Standard output/error logs, data files

## Interfaces & Ecosystem
- **ASE**: Often used in conjunction with ASE scripts
- **GPAW**: Commonly used by the GPAW community

## Workflow and Usage
1. Create folder structure: `simulations/run1/`
2. Place script `run.py` in folder.
3. Run `mq submit simulations/run1/`
4. Monitor with `mq list`.
5. If failed, fix and `mq resubmit`.

## Performance Characteristics
- Very low overhead
- Simple file-system based state
- Ideal for personal high-throughput management

## Application Areas
- Managing batches of DFT calculations
- Parameter sweeps
- Simple dependency chains (relax -> band structure)

## Community and Support
- Developed at DTU Physics (Jakob Schiøtz and contributors)
- Used by CAMD/Thygesen groups

## Verification & Sources
**Primary sources**:
1. Documentation: https://myqueue.readthedocs.io/
2. GitLab: https://gitlab.com/schmidt-group/myqueue

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitLab)
- Development: ACTIVE
- Applications: Simple job management, folder-based workflows
