# CSPML

## Overview
CSPML (Crystal Structure Prediction with Machine Learning) is a template-based crystal structure prediction method using metric learning for element substitution. It predicts stable structures by selecting templates from known crystal structures.

## Theoretical Basis
- Template-based structure prediction
- Metric learning for element substitution
- Crystal structure similarity
- Formation energy prediction
- Template selection optimization

## Key Capabilities
- Element substitution-based CSP
- Metric learning for template selection
- No ab initio calculations needed
- Fast structure prediction
- Wide chemical space coverage

**Sources**: Comp. Mater. Sci. 211, 111496 (2022)

## Key Strengths

### Speed:
- No DFT required
- Fast prediction
- Large-scale screening

### Methodology:
- Metric learning
- Template selection
- Element substitution

### Coverage:
- Wide chemical space
- Multiple crystal systems
- General applicability

## Inputs & Outputs
- **Input formats**: Chemical composition
- **Output data types**: Predicted structures, template matches

## Interfaces & Ecosystem
- **Python**: TensorFlow-based
- **Databases**: Template database
- **Validation**: DFT codes

## Workflow and Usage
1. Input chemical composition
2. Run metric learning model
3. Select best templates
4. Generate substituted structures
5. Validate with DFT (optional)

## Performance Characteristics
- Very fast prediction
- 50-65% accuracy on crystal systems
- Scales to large datasets

## Computational Cost
- Minimal (ML inference)
- No DFT required
- Fast screening

## Best Practices
- Use diverse template database
- Validate top predictions
- Consider multiple templates
- Check structural validity

## Limitations & Known Constraints
- Template-dependent
- May miss novel structures
- Accuracy varies by system

## Application Areas
- High-throughput screening
- Initial structure guessing
- Materials discovery
- Database expansion

## Comparison with Other Codes
- **vs TCSP**: Similar approach, different ML
- **vs USPEX**: CSPML faster, less accurate
- **Unique strength**: Fast template-based prediction, no DFT

## Community and Support
- Open-source (GitHub)
- Published methodology
- Code available

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/Minoru938/CSPML
2. Publication: Comp. Mater. Sci. 211, 111496 (2022)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: MAINTAINED
- Applications: Fast CSP, template-based prediction
