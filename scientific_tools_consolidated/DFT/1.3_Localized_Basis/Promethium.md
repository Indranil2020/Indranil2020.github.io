# Promethium

## Official Resources
- Homepage: https://qcware.com/promethium
- Documentation: Commercial / QC Ware
- Source Repository: Proprietary (SaaS)
- License: Commercial (AWS Marketplace)

## Overview
Promethium is a cloud-native, GPU-accelerated quantum chemistry platform developed by QC Ware. It is delivered as a Software-as-a-Service (SaaS), primarily via AWS. Promethium utilizes advanced algorithms optimized for NVIDIA GPUs (H100/A100) to perform Density Functional Theory (DFT) calculations with exceptional speed and throughput, handling systems up to 2,000 atoms.

**Scientific domain**: Quantum Chemistry, High-Throughput Screening, DFT  
**Target user community**: Pharma, Materials Enterprise, R&D

## Theoretical Methods
- Density Functional Theory (DFT)
- Hybrid Functionals (B3LYP, PBE0, etc.)
- Basis sets (Gaussian type)
- Geometry Optimization
- Vibrational Frequency analysis
- Solvation models
- Transition State search

## Capabilities (CRITICAL)
- GPU-accelerated DFT
- High-throughput batch processing
- Large system support (2000 atoms)
- Conformer generation and ranking
- Reaction path analysis
- Cloud-native deployment
- Seamless scaling
- Python API integration

## Key Strengths

### Speed:
- 50x-100x faster than traditional CPU codes
- GPU-native algorithms
- Fast hybrid functionals
- Rapid turnaround for large batches

### Scale:
- Handles protein-ligand pockets
- Large nanoclusters
- Supramolecular complexes
- Massive datasets

### Ease of Use:
- Cloud managed (no cluster maintenance)
- Python SDK
- Integration with workflows
- On-demand resources

## Inputs & Outputs
- **Input**: Molecular structures, workflow configuration
- **Output**: Energies, Properties, Geometries
- **API**: Modern Python interface

## Interfaces & Ecosystem
- **AWS**: Deployed on Amazon Web Services
- **NVIDIA**: Optimized for Tensor Cores
- **Python**: Client library for submission/retrieval

## Advanced Features

### GPU Algorithms:
- Fast exchange-correlation evaluation
- Optimized grid integration
- Direct SCF schemes
- Memory management for large systems

### Workflows:
- Automated conformer search
- Reaction profile scans
- High-throughput screening pipelines

## Performance Characteristics
- **Speed**: State-of-the-art GPU DFT
- **Accuracy**: Standard DFT precision
- **System size**: 100-2000 atoms
- **Scaling**: Linear/near-linear for many tasks

## Computational Cost
- **Model**: Pay-per-compute / Subscription
- **Efficiency**: Highly efficient per calculation
- **Overhead**: Cloud latency (minimal)

## Limitations & Known Constraints
- **Access**: Commercial SaaS only
- **Cost**: Cloud costs apply
- **Customization**: Less flexible than open source code
- **Methods**: Focused on standard DFT

## Comparison with Other Codes
- **vs Gaussian/ORCA**: Promethium is cloud-native GPU SaaS
- **vs TeraChem**: Similar GPU focus, Promethium is SaaS model
- **vs CPU codes**: Significantly faster for large hybrid DFT
- **Unique strength**: Turn-key GPU DFT in the cloud

## Application Areas
### Computational Chemistry R&D:
- **Binding Free Energies**: Rapid calculation of GPCR-ligand interactions
- **Transition State Finding**: Automated location of catalytic barriers
- **Conformational Search**: Exhaustive sampling of flexible drug-like molecules
- **Reaction Networks**: Mapping complex organic reaction pathways

### Materials Engineering:
- **Battery Materials**: Screening electrolyte decomposition pathways
- **OLED Design**: Excited state properties (TD-DFT) for emitters (future)
- **Catalysis**: High-throughput screening of organometallic catalysts
- **Nanomaterials**: Structure prediction of large metal nanoclusters

## Best Practices
### Workflow Optimization:
- **Batch Submission**: Group calculations into large batches to amortize network latency
- **Resource Selection**: Choose appropriate instance types (H100 vs A100) based on system size
- **Storage Management**: Periodically archive or download results to avoid storage costs
- **Monitoring**: Use dashboard to track compute usage and budget

### Calculation Settings:
- **Basis Sets**: Use def2-SVP/TZVP for optimal GPU performance
- **Grids**: Standard grids (SG-1/SG-2) are highly optimized
- **Convergence**: Enable level-shifting for difficult metalloproteins
- **Geometry**: Pre-optimize with clean force fields (MMFF94) before DFT

## Community and Support
- QC Ware support team
- AWS Marketplace support

## Verification & Sources
**Primary sources**:
1. https://qcware.com/promethium
2. Press releases (NVIDIA, AWS)
3. QC Ware publications

**Confidence**: VERIFIED
- Status: Active commercial product
- Tech: Verified GPU acceleration
