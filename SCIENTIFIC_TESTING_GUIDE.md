# Scientific Testing Framework Guide

This guide explains the scientific testing framework implemented in Phase 4 of the CryoProtect project. The framework provides a robust way to verify the scientific accuracy of the various models in the system.

## Overview

The scientific testing framework is designed to compare the predictions of scientific models against known reference values from literature or experimental data. This ensures that the models are scientifically accurate and reliable.

The framework includes:

1. A set of test cases for different model types
2. Reference data with known values
3. Tolerance-based comparison methods
4. A structured approach to validating scientific accuracy

## Test Structure

The test suite is organized by model type:

1. **ConcentrationModelAccuracyTestCase**: Tests for concentration-dependent models
2. **TemperatureModelAccuracyTestCase**: Tests for temperature-dependent models
3. **GlassTransitionAccuracyTestCase**: Tests for glass transition temperature prediction
4. **SynergyPredictionAccuracyTestCase**: Tests for synergy prediction between components

Each test case loads reference data and compares model predictions against expected values within specified tolerances.

## Reference Data

The reference data is stored in JSON format in `tests/reference_data/scientific_reference_data.json`. This data includes:

- **Concentration Models**: Test points with concentrations and expected property values
- **Temperature Models**: Test points with temperatures and expected property values
- **Glass Transition**: Expected glass transition temperatures for pure compounds and mixtures
- **Synergy Prediction**: Expected synergy score ranges for different component combinations

Example reference data:

```json
{
  "concentration_models": {
    "linear": [
      {
        "molecule_id": "glycerol",
        "parameters": {
          "slope": 0.5,
          "intercept": 1.0
        },
        "test_points": [
          {"concentration": 0.0, "expected_value": 1.0},
          {"concentration": 1.0, "expected_value": 1.5}
        ]
      }
    ]
  }
}
```

## Testing Methodology

The framework employs the following methodology:

1. **Load Reference Data**: Test cases load reference data from JSON files
2. **Create Model Instances**: Models are created with parameters from reference data
3. **Run Model Calculations**: Models calculate predictions for test points
4. **Compare Results**: Predictions are compared to expected values with tolerances
5. **Assert Accuracy**: Tests pass if predictions are within tolerance ranges

### Tolerance-Based Comparisons

The framework uses a tolerance-based approach to account for numerical differences and implementation variations:

```python
def assert_close(self, expected: float, actual: float, 
                 tolerance: float = TOLERANCE, msg: str = None) -> None:
    # Calculate absolute tolerance based on expected value
    abs_tolerance = abs(expected * tolerance)
    
    # Use built-in assertAlmostEqual with absolute tolerance
    self.assertAlmostEqual(
        expected, actual, delta=abs_tolerance,
        msg=msg or f"Expected {expected} but got {actual} (tolerance: {tolerance*100}%)"
    )
```

By default, a 5% tolerance is used, but this can be adjusted for specific tests as needed.

## Running the Tests

To run the scientific accuracy tests:

```bash
# Run all scientific accuracy tests
python -m unittest tests/test_scientific_accuracy.py

# Run a specific test case
python -m unittest tests.test_scientific_accuracy.ConcentrationModelAccuracyTestCase

# Run a specific test method
python -m unittest tests.test_scientific_accuracy.ConcentrationModelAccuracyTestCase.test_linear_concentration_model
```

## Extending the Framework

### Adding New Reference Data

To add new reference data:

1. **Identify Reliable Sources**: Use peer-reviewed literature or experimental data
2. **Format Data**: Structure the data according to the existing JSON format
3. **Add to Reference File**: Add the data to `scientific_reference_data.json` or create a new file

### Adding New Test Cases

To add new test cases:

1. **Create Test Class**: Extend `ScientificAccuracyTestCase` for the new model type
2. **Implement Test Methods**: Add test methods for different aspects of the model
3. **Register Reference Data**: Ensure reference data is available for the tests
4. **Validate**: Run the tests to ensure they correctly validate the model

## Maintenance

To maintain the testing framework:

1. **Update Reference Data**: As new scientific literature becomes available
2. **Review Tolerances**: Adjust tolerances based on model refinements
3. **Add New Model Tests**: When new scientific models are implemented
4. **Document Changes**: Keep documentation current with framework evolution

## Best Practices

1. **Use Reliable Sources**: Only use peer-reviewed or experimentally validated reference data
2. **Test Edge Cases**: Include extreme values to test model robustness
3. **Validate Implementation**: Verify that models correctly implement underlying equations
4. **Document Limitations**: Note known limitations or assumptions in tests
5. **Regular Testing**: Run tests after any model changes

## Conclusion

The scientific testing framework provides a robust foundation for validating the accuracy and reliability of the scientific models in the CryoProtect system. By comparing model predictions against known reference values, we can ensure that our models correctly implement scientific principles and provide accurate predictions for cryoprotectant behavior.