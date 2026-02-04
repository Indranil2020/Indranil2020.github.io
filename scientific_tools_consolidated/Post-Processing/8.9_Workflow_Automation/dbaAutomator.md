# dbaAutomator

## Official Resources
- **Homepage**: https://github.com/xingyu-alfred-liu/dbaAutomator
- **BerkeleyGW News**: https://berkeleygw.org/2020/01/29/dbaautomator-now-available-for-berkeleygw/
- **Source Repository**: https://github.com/xingyu-alfred-liu/dbaAutomator
- **License**: MIT License (implied by GitHub usage, check repo)

## Overview
dbaAutomator (Double-Bader Analysis Automator) is a Python-based tool designed to assist users of the **BerkeleyGW** package. It automates the workflow for performing **Double-Bader Analysis (DBA)** on excitons in molecular crystals. This analysis characterizes the hole and electron distributions of excitons, quantifying the degree of charge transfer. It also provides functionality to verify the convergence of fine k-point grids for GW-BSE calculations.

**Scientific domain**: Exciton analysis, Charge transfer, GW-BSE post-processing
**Target user community**: BerkeleyGW users, researchers studying organic photovoltaics and molecular crystals

## Capabilities
- **Double-Bader Analysis**: Automates the integration of exciton probability densities over Bader volumes to classify excitons (Frenkel, Charge Transfer, etc.).
- **Convergence Testing**: Verifies convergence of the fine k-point grid required for BSE calculations.
- **Automation**: Handles the generation of input files and parsing of outputs for the Bader analysis steps.

## Inputs & Outputs
- **Inputs**: BerkeleyGW output files (exciton wavefunctions), Bader charge analysis outputs.
- **Outputs**: Charge transfer characterization (CT numbers), convergence plots.

## Interfaces & Ecosystem
- **BerkeleyGW**: The primary code for generating the excitonic data.
- **Bader**: Uses the Henkelman Group's Bader charge analysis code for the volume integration.
- **Python**: Written in Python.

## Verification & Sources
- **Primary Source**: [GitHub Repository](https://github.com/xingyu-alfred-liu/dbaAutomator)
- **Announcement**: [BerkeleyGW News](https://berkeleygw.org/2020/01/29/dbaautomator-now-available-for-berkeleygw/)
- **Publication**: X. Liu et al., *Journal of Physics: Condensed Matter* 32, 185901 (2020). [DOI: 10.1088/1361-648X/ab699e](https://iopscience.iop.org/article/10.1088/1361-648X/ab699e)
