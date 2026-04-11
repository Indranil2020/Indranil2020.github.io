# AiiDA Plugin Registry

## Official Resources
- Homepage: https://aiidateam.github.io/aiida-registry/
- Documentation: https://aiidateam.github.io/aiida-registry/
- Source Repository: https://github.com/aiidateam/aiida-registry
- License: MIT License

## Overview
The AiiDA Plugin Registry is not a software tool itself but the central catalog of all available plugins for the AiiDA ecosystem. It lists verified and community-contributed plugins that interface AiiDA with various simulation codes, schedulers, and data tools. It serves as the primary discovery mechanism for AiiDA users to find extensions for their specific codes.

**Scientific domain**: Software catalog, ecosystem management  
**Target user community**: AiiDA users looking for code support

## Capabilities (CRITICAL)
- **Discovery**: Searchable list of plugins (e.g., VASP, QE, CP2K, LAMMPS).
- **Metadata**: Provides installation info, documentation links, and PyPI package names.
- **Status**: Indicates development status (stable, beta, development) and CI status.
- **Registration**: Process for developers to register new plugins.

**Sources**: AiiDA Plugin Registry website

## Inputs & Outputs
- **Input**: Pull requests to register new plugins
- **Output**: JSON database of plugins, web interface

## Interfaces & Ecosystem
- **AiiDA**: The core framework
- **PyPI**: Where plugins are hosted
- **GitHub**: Where source code lives

## Workflow and Usage
1. Visit https://aiidateam.github.io/aiida-registry/
2. Search for a code (e.g., "SIESTA").
3. Find the plugin package name (e.g., `aiida-siesta`).
4. Install: `pip install aiida-siesta`.

## Performance Characteristics
- Static website / JSON API
- Automatically updated via GitHub Actions

## Application Areas
- Finding tools for AiiDA workflows
- Ecosystem health monitoring

## Community and Support
- Maintained by AiiDA team
- Community contributions for new plugins

## Verification & Sources
**Primary sources**:
1. Homepage: https://aiidateam.github.io/aiida-registry/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Plugin discovery
