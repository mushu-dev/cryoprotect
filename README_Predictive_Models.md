# CryoProtect Analyzer - Predictive Models

This document provides an overview of the predictive models functionality in the CryoProtect Analyzer application.

## Overview

The predictive models module allows users to train machine learning models to predict cryoprotectant performance based on molecular properties. These models can be used to predict the effectiveness of new cryoprotectant mixtures without the need for experimental testing, saving time and resources.

## Features

- **Multiple Prediction Algorithms**: The system supports various machine learning algorithms, including:
  - Linear Regression
  - Ridge Regression
  - Lasso Regression
  - Random Forest
  - Gradient Boosting
  - Neural Networks

- **Model Training**: Train models using experimental data stored in the system. The models learn to predict cryoprotectant properties based on molecular features.

- **Confidence Intervals**: All predictions include confidence intervals to indicate the reliability of the prediction.

- **Feature Importance**: Understand which molecular properties have the most significant impact on cryoprotection effectiveness.

- **Model Management**: View, compare, and delete trained models.

## How It Works

### Molecular Features

The predictive models use the following molecular features:

1. **Hydrogen Bonding**: Number of hydrogen bond donors and acceptors
2. **LogP**: Partition coefficient (measure of hydrophobicity/hydrophilicity)
3. **TPSA**: Topological polar surface area
4. **Molecular Properties**: Weight, size, rotatable bonds, etc.
5. **Functional Groups**: Presence and count of important functional groups (hydroxyl, ether, etc.)
6. **Permeability**: Estimated cell membrane permeability

### Training Process

1. The system collects experimental data from the database
2. Features are extracted from the molecular structures
3. The model is trained using the selected algorithm
4. Model performance is evaluated using validation metrics
5. The trained model is saved for future use

### Prediction Process

1. For a given mixture, the system extracts features from each component
2. The model predicts the property value for each component
3. Component predictions are combined using concentration-weighted averaging
4. A synergy bonus is applied for mixtures with complementary properties
5. The final prediction is returned with confidence intervals

## Using the Predictive Models

### Training a Model

1. Navigate to the "Predictive Models" page
2. Select the "Train Model" tab
3. Enter the property name (e.g., "Cryoprotection Score")
4. Select an algorithm
5. Configure hyperparameters (optional)
6. Click "Train Model"

### Making Predictions

1. Navigate to the "Predictive Models" page
2. Select the "Make Prediction" tab
3. Select a mixture
4. Select an algorithm
5. Click "Make Prediction"
6. View the prediction results, including confidence intervals

## Integration with Other Modules

The predictive models functionality integrates with other modules in the CryoProtect Analyzer:

- **Molecules Module**: Uses molecular properties calculated by RDKit
- **Mixtures Module**: Makes predictions for mixtures of cryoprotectants
- **Experiments Module**: Uses experimental data for model training
- **Comparisons Module**: Compares predictions with experimental results

## Technical Implementation

The predictive models functionality is implemented in the following files:

- `api/predictive_models.py`: Core functionality for model training and prediction
- `api/predictive_models_resources.py`: API endpoints for the predictive models
- `static/js/predictive-models.js`: Frontend JavaScript for the predictive models page
- `templates/predictive_models.html`: HTML template for the predictive models page
- `tests/test_predictive_models.py`: Unit tests for the predictive models functionality

## Best Practices

- **Data Quality**: Ensure experimental data is accurate and consistent
- **Feature Selection**: Consider which molecular properties are most relevant for your prediction task
- **Algorithm Selection**: Different algorithms may perform better for different properties
- **Hyperparameter Tuning**: Optimize model hyperparameters for better performance
- **Validation**: Always check prediction confidence and compare with experimental results when available

## Future Enhancements

Potential future enhancements to the predictive models functionality:

1. **Advanced Algorithms**: Support for more advanced machine learning algorithms
2. **Automated Hyperparameter Optimization**: Automatic tuning of model hyperparameters
3. **Ensemble Models**: Combining multiple models for improved prediction accuracy
4. **Transfer Learning**: Using pre-trained models for related properties
5. **Explainable AI**: Better visualization and explanation of model predictions
6. **Active Learning**: Suggesting experiments to improve model performance
7. **Molecular Fingerprints**: Using molecular fingerprints as additional features
8. **3D Structure Analysis**: Incorporating 3D structural information for better predictions

---

## API Endpoint: Predictive Scoring for Cryoprotectant Effectiveness

### Endpoint

```
GET /mixtures/<mixture_id>/predict
```

- **Description:** Returns a predictive effectiveness score for a cryoprotectant mixture, including confidence and contributing factors.
- **Query Parameters:**
  - `algorithm` (optional): The prediction algorithm to use (e.g., `random_forest`, `ridge_regression`). Default: `random_forest`.

### Example Request

```
GET /mixtures/12345/predict?algorithm=random_forest
```

### Example Response

```json
{
  "property_name": "Cryoprotection Score",
  "prediction": 82.5,
  "confidence": 5.2,
  "confidence_interval": [77.3, 87.7],
  "algorithm": "random_forest",
  "algorithm_name": "Random Forest",
  "feature_importance": {
    "logp": 0.18,
    "tpsa": 0.12,
    "h_bond_donors": 0.09,
    "molecular_weight": 0.15,
    "...": "..."
  }
}
```

- `prediction`: The estimated effectiveness score (higher is better).
- `confidence`: The model's estimated uncertainty (lower is better).
- `confidence_interval`: The range in which the true value is likely to fall.
- `feature_importance`: Key molecular features and their relative influence on the prediction.

### Interpretation

- **Effectiveness Score:** A higher score indicates a more effective cryoprotectant or mixture, as predicted by the model.
- **Confidence:** Lower values indicate higher certainty in the prediction.
- **Feature Importance:** Shows which molecular properties most influenced the prediction, aiding scientific interpretation and transparency.

---

## Scientific Basis and Assumptions

- **Data-Driven:** The predictive scoring algorithm is trained on experimental data linking molecular/mixture features to cryoprotectant effectiveness.
- **Feature Engineering:** Uses molecular descriptors (hydrogen bonding, LogP, TPSA, etc.) and mixture composition.
- **Mixture Scoring:** Predictions for mixtures are calculated as a concentration-weighted average of component predictions, with a synergy bonus for mixtures with complementary properties.
- **Model Selection:** Multiple machine learning algorithms are supported; Random Forest is the default for robustness.
- **Interpretability:** Feature importances are provided to explain which properties most influenced the score.
- **Assumptions:** The model assumes that the training data is representative and that the selected features capture the main determinants of cryoprotectant effectiveness. Synergy effects are modeled simply and may not capture all real-world interactions.

---

## Frontend Integration

- The `/mixtures/<mixture_id>/predict` endpoint can be called from the frontend to display predictive scores, confidence, and contributing factors.
- See `static/js/predictive-models.js` for usage examples.

---

## Contact

For questions or suggestions regarding the predictive scoring algorithm, please contact the CryoProtect Analyzer development team.