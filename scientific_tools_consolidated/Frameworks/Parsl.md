# Parsl (Parallel Scripting Library)

## Official Resources
- Homepage: https://parsl-project.org/
- Documentation: https://parsl.readthedocs.io/
- Source Repository: https://github.com/Parsl/parsl
- License: Apache License 2.0

## Overview
Parsl is a Python library for flexible, parallel scripting. It allows researchers to build parallel applications composed of Python functions and external components (executables) that can run on arbitrary computing resources, from laptops to supercomputers. It abstracts the execution model, allowing the same script to scale from a single core to thousands of nodes.

**Scientific domain**: Parallel computing, workflow management, HPC scripting  
**Target user community**: Researchers needing scalable parallelism in Python

## Capabilities (CRITICAL)
- **Elastic Parallelism**: Dynamically scales resources based on workload.
- **DataFlow Kernel**: Asynchronous execution based on data dependencies (Futures).
- **Providers**: Interfaces for SLURM, PBS, Cobalt, Kubernetes, AWS, etc.
- **Executors**: High-performance execution engines (HighThroughputExecutor).
- **Pythonic**: Decorator-based syntax (`@python_app`, `@bash_app`).
- **Monitoring**: Real-time monitoring database.

**Sources**: Parsl website, HPDC '19

## Inputs & Outputs
- **Input formats**: Python scripts
- **Output data types**: Python objects (Futures), files

## Interfaces & Ecosystem
- **Jupyter**: Works well within notebooks
- **HPC**: Deep integration with supercomputing schedulers (Argonne, NERSC, etc.)
- **FuncX**: Related project for function-as-a-service

## Workflow and Usage
1. Configure executor (e.g., local threads or SLURM config).
2. Decorate functions:
   ```python
   @python_app
   def simulate(x):
       return x**2
   ```
3. Call functions (returns Future).
4. `result = simulate(10).result()` (blocks until done).

## Performance Characteristics
- High throughput (thousands of tasks/sec with HighThroughputExecutor)
- Low latency
- Scalable to extreme-scale systems (Exascale)

## Application Areas
- High-energy physics analyses
- Cosmology simulations
- Materials science screening
- Bioinformatics

## Community and Support
- Developed by University of Chicago / Argonne National Laboratory
- Active academic project
- Annual workshops

## Verification & Sources
**Primary sources**:
1. Homepage: https://parsl-project.org/
2. GitHub: https://github.com/Parsl/parsl

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Parallel scripting, HPC workflows
