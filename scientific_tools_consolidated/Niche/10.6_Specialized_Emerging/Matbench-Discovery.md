# Matbench-Discovery

## Official Resources
- Source Repository: https://github.com/janosh/matbench-discovery
- Website: https://matbench-discovery.materialsproject.org/
- License: Open source (MIT)

## Overview
**Matbench-Discovery** is an interactive leaderboard and evaluation framework for ML models simulating high-throughput materials discovery. It ranks 20+ models on stability prediction, structure relaxation, and thermal conductivity tasks.

**Scientific domain**: Benchmarking framework for MLIP and materials discovery models  
**Target user community**: Researchers evaluating and comparing MLIP models

## Theoretical Methods
- Multi-task benchmarking
- Stability prediction (F1, Precision, Recall)
- Structure relaxation (energy above hull)
- Thermal conductivity prediction
- Model comparison metrics

## Capabilities (CRITICAL)
- 20+ model leaderboard
- Stability prediction benchmark
- Structure relaxation benchmark
- Thermal conductivity benchmark
- Interactive website
- Reproducible evaluation

**Sources**: GitHub repository

## Key Strengths

### Benchmarking:
- Fair model comparison
- Multiple tasks
- Standardized evaluation
- Reproducible results

### Leaderboard:
- Interactive website
- Real-time updates
- Model metadata
- Community contributions

### Comprehensive:
- 20+ models ranked
- MACE, CHGNet, M3GNet, SevenNet, etc.
- Multiple metrics
- Statistical analysis

## Inputs & Outputs
- **Input formats**: Model predictions
- **Output data types**: Rankings, metrics, analysis

## Interfaces & Ecosystem
- **Python**: Core
- **Website**: Interactive
- **Materials Project**: Data source

## Performance Characteristics
- **Speed**: Evaluation takes minutes
- **Accuracy**: N/A (benchmarking tool)
- **System size**: N/A
- **Automation**: Full

## Computational Cost
- **Evaluation**: Minutes per model
- **No training**: Uses existing predictions

## Limitations & Known Constraints
- **PBE reference**: Benchmarked against PBE data
- **Limited tasks**: 3 main tasks
- **Model-dependent**: Results depend on model quality
- **Computational cost**: Running models is expensive

## Comparison with Other Codes
- **vs Matbench**: Matbench-Discovery is UIP-focused, Matbench is general
- **vs MLIP Arena**: Matbench-Discovery is materials, MLIP Arena is general
- **Unique strength**: Interactive leaderboard ranking 20+ UIP models on materials discovery tasks

## Application Areas

### Model Selection:
- Choose best UIP for application
- Compare model accuracy
- Track model improvements
- Guide development

### Research:
- Benchmark new models
- Analyze model strengths
- Identify failure modes

## Best Practices
- Check leaderboard before choosing model
- Compare multiple models
- Consider task-specific performance
- Reproduce results

## Community and Support
- Open source (MIT)
- Materials Project affiliated
- Interactive website
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/janosh/matbench-discovery

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Interactive leaderboard ranking 20+ UIP models on materials discovery tasks
