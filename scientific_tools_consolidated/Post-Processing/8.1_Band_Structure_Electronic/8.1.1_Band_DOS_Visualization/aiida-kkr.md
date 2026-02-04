# aiida-kkr

## Official Resources
- **GitHub**: https://github.com/JuDFTteam/aiida-kkr
- **Documentation**: https://aiida-kkr.readthedocs.io/
- **License**: MIT License

## Overview
aiida-kkr is an AiiDA plugin for running KKR (Korringa-Kohn-Rostoker) Green's function calculations, providing workflow automation for KKR-based electronic structure calculations.

**Scientific domain**: KKR Green's function, workflow automation
**Target user community**: KKR users needing automated workflows

## Capabilities (CRITICAL)
- **KKR Calculations**: Green's function method
- **Workflow Automation**: AiiDA integration
- **Impurity Calculations**: Defect modeling
- **Transport**: Ballistic transport

## Key Strengths
- AiiDA workflow management
- JuKKR code interface
- Impurity embedding
- High-throughput support

## Requirements
- AiiDA framework
- JuKKR code

## Installation
```bash
pip install aiida-kkr
```

## Limitations & Known Constraints
- **KKR-specific**: Only interfaces with JuKKR codes
- **AiiDA dependency**: Requires AiiDA framework setup
- **Learning curve**: AiiDA workflow and KKR method knowledge needed
- **Specialized**: Green's function method less common than plane-wave

## Comparison with Other Tools
- **vs aiida-fleur**: aiida-kkr for KKR, aiida-fleur for FLEUR
- **vs masci-tools**: aiida-kkr workflows, masci-tools post-processing
- **Unique strength**: AiiDA workflow automation for KKR Green's function calculations

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Developer: JuDFTteam
