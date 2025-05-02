# CryoProtect v2 API Integration Documentation

## Overview

This document provides comprehensive documentation for the CryoProtect v2 API, which integrates with the Supabase database. The API provides endpoints for managing molecules, mixtures, predictions, experiments, and more.

## Base URL

All API endpoints are relative to the base URL: `http://localhost:5000/api/v1`

## Authentication

Most endpoints require authentication using a JWT token. Include the token in the Authorization header:

```
Authorization: Bearer <your_token>
```

To obtain a token, use the `/auth/login` endpoint.

## Error Handling

All endpoints return standard HTTP status codes:

- 200: Success
- 201: Created
- 400: Bad Request
- 401: Unauthorized
- 403: Forbidden
- 404: Not Found
- 409: Conflict
- 500: Internal Server Error

Error responses include a JSON object with a message field:

```json
{
  "message": "Error message details"
}
```

## Endpoints

### Molecules

#### List Molecules

```
GET /molecules
```

Query Parameters:
- `limit` (optional): Maximum number of molecules to return (default: 100, max: 500)
- `offset` (optional): Number of molecules to skip (default: 0)

Response:
```json
[
  {
    "id": "uuid",
    "cid": 123456,
    "name": "Ethanol",
    "molecular_formula": "C2H6O",
    "smiles": "CCO",
    "pubchem_link": "https://pubchem.ncbi.nlm.nih.gov/compound/123456",
    "created_at": "2025-04-18T15:00:00Z",
    "updated_at": "2025-04-18T15:00:00Z",
    "properties": {
      "molecular_weight": 46.07,
      "logp": -0.31,
      "tpsa": 20.23
    }
  }
]
```

#### Get Molecule

```
GET /molecules/{molecule_id}
```

Path Parameters:
- `molecule_id`: UUID of the molecule

Response:
```json
{
  "id": "uuid",
  "cid": 123456,
  "name": "Ethanol",
  "molecular_formula": "C2H6O",
  "smiles": "CCO",
  "pubchem_link": "https://pubchem.ncbi.nlm.nih.gov/compound/123456",
  "created_at": "2025-04-18T15:00:00Z",
  "updated_at": "2025-04-18T15:00:00Z",
  "properties": {
    "molecular_weight": 46.07,
    "logp": -0.31,
    "tpsa": 20.23
  }
}
```

#### Create Molecule

```
POST /molecules
```

Request Body:
```json
{
  "cid": 123456
}
```

Response:
```json
{
  "id": "uuid",
  "cid": 123456,
  "name": "Ethanol",
  "molecular_formula": "C2H6O",
  "smiles": "CCO",
  "pubchem_link": "https://pubchem.ncbi.nlm.nih.gov/compound/123456",
  "created_at": "2025-04-18T15:00:00Z",
  "updated_at": "2025-04-18T15:00:00Z",
  "properties": {
    "molecular_weight": 46.07,
    "logp": -0.31,
    "tpsa": 20.23
  }
}
```

#### Update Molecule

```
PUT /molecules/{molecule_id}
```

Path Parameters:
- `molecule_id`: UUID of the molecule

Request Body:
```json
{
  "name": "Updated Ethanol",
  "formula": "C2H6O",
  "smiles": "CCO"
}
```

Response:
```json
{
  "id": "uuid",
  "cid": 123456,
  "name": "Updated Ethanol",
  "molecular_formula": "C2H6O",
  "smiles": "CCO",
  "pubchem_link": "https://pubchem.ncbi.nlm.nih.gov/compound/123456",
  "created_at": "2025-04-18T15:00:00Z",
  "updated_at": "2025-04-18T15:30:00Z",
  "properties": {
    "molecular_weight": 46.07,
    "logp": -0.31,
    "tpsa": 20.23
  }
}
```

#### Delete Molecule

```
DELETE /molecules/{molecule_id}
```

Path Parameters:
- `molecule_id`: UUID of the molecule

Response:
```json
{
  "message": "Molecule with ID {molecule_id} deleted successfully"
}
```

### Mixtures

#### List Mixtures

```
GET /mixtures
```

Query Parameters:
- `limit` (optional): Maximum number of mixtures to return (default: 100, max: 500)
- `offset` (optional): Number of mixtures to skip (default: 0)

Response:
```json
[
  {
    "id": "uuid",
    "name": "Ethanol-Glycerol Mixture",
    "description": "A mixture of ethanol and glycerol",
    "created_at": "2025-04-18T15:00:00Z",
    "updated_at": "2025-04-18T15:00:00Z",
    "components": [
      {
        "molecule_id": "uuid",
        "name": "Ethanol",
        "concentration": 70,
        "concentration_unit": "%v/v"
      },
      {
        "molecule_id": "uuid",
        "name": "Glycerol",
        "concentration": 30,
        "concentration_unit": "%v/v"
      }
    ]
  }
]
```

#### Get Mixture

```
GET /mixtures/{mixture_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture

Response:
```json
{
  "id": "uuid",
  "name": "Ethanol-Glycerol Mixture",
  "description": "A mixture of ethanol and glycerol",
  "created_at": "2025-04-18T15:00:00Z",
  "updated_at": "2025-04-18T15:00:00Z",
  "components": [
    {
      "molecule_id": "uuid",
      "name": "Ethanol",
      "concentration": 70,
      "concentration_unit": "%v/v"
    },
    {
      "molecule_id": "uuid",
      "name": "Glycerol",
      "concentration": 30,
      "concentration_unit": "%v/v"
    }
  ]
}
```

#### Create Mixture

```
POST /mixtures
```

Request Body:
```json
{
  "name": "Ethanol-Glycerol Mixture",
  "description": "A mixture of ethanol and glycerol",
  "components": [
    {
      "molecule_id": "uuid",
      "concentration": 70,
      "concentration_unit": "%v/v"
    },
    {
      "molecule_id": "uuid",
      "concentration": 30,
      "concentration_unit": "%v/v"
    }
  ]
}
```

Response:
```json
{
  "id": "uuid",
  "name": "Ethanol-Glycerol Mixture",
  "description": "A mixture of ethanol and glycerol",
  "created_at": "2025-04-18T15:00:00Z",
  "updated_at": "2025-04-18T15:00:00Z",
  "components": [
    {
      "molecule_id": "uuid",
      "name": "Ethanol",
      "concentration": 70,
      "concentration_unit": "%v/v"
    },
    {
      "molecule_id": "uuid",
      "name": "Glycerol",
      "concentration": 30,
      "concentration_unit": "%v/v"
    }
  ]
}
```

#### Update Mixture

```
PUT /mixtures/{mixture_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture

Request Body:
```json
{
  "name": "Updated Ethanol-Glycerol Mixture",
  "description": "An updated mixture of ethanol and glycerol",
  "components": [
    {
      "molecule_id": "uuid",
      "concentration": 60,
      "concentration_unit": "%v/v"
    },
    {
      "molecule_id": "uuid",
      "concentration": 40,
      "concentration_unit": "%v/v"
    }
  ]
}
```

Response:
```json
{
  "id": "uuid",
  "name": "Updated Ethanol-Glycerol Mixture",
  "description": "An updated mixture of ethanol and glycerol",
  "created_at": "2025-04-18T15:00:00Z",
  "updated_at": "2025-04-18T15:30:00Z",
  "components": [
    {
      "molecule_id": "uuid",
      "name": "Ethanol",
      "concentration": 60,
      "concentration_unit": "%v/v"
    },
    {
      "molecule_id": "uuid",
      "name": "Glycerol",
      "concentration": 40,
      "concentration_unit": "%v/v"
    }
  ]
}
```

#### Delete Mixture

```
DELETE /mixtures/{mixture_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture

Response:
```json
{
  "message": "Mixture with ID {mixture_id} deleted successfully"
}
```

### Predictions

#### List Predictions for a Mixture

```
GET /mixtures/{mixture_id}/predictions
```

Path Parameters:
- `mixture_id`: UUID of the mixture

Query Parameters:
- `limit` (optional): Maximum number of predictions to return (default: 100, max: 500)
- `offset` (optional): Number of predictions to skip (default: 0)

Response:
```json
[
  {
    "id": "uuid",
    "mixture_id": "uuid",
    "property_type_id": "uuid",
    "property_name": "Glass Transition Temperature",
    "numeric_value": -45.2,
    "text_value": null,
    "boolean_value": null,
    "confidence": 0.85,
    "calculation_method": "ML Model v1",
    "created_at": "2025-04-18T15:00:00Z"
  }
]
```

#### Get Prediction

```
GET /mixtures/{mixture_id}/predictions/{prediction_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture
- `prediction_id`: UUID of the prediction

Response:
```json
{
  "id": "uuid",
  "mixture_id": "uuid",
  "property_type_id": "uuid",
  "property_name": "Glass Transition Temperature",
  "numeric_value": -45.2,
  "text_value": null,
  "boolean_value": null,
  "confidence": 0.85,
  "calculation_method": "ML Model v1",
  "created_at": "2025-04-18T15:00:00Z"
}
```

#### Create Prediction

```
POST /mixtures/{mixture_id}/predictions
```

Path Parameters:
- `mixture_id`: UUID of the mixture

Request Body:
```json
{
  "property_name": "Glass Transition Temperature",
  "value": -45.2,
  "confidence": 0.85,
  "calculation_method": "ML Model v1"
}
```

Response:
```json
{
  "id": "uuid",
  "mixture_id": "uuid",
  "property_type_id": "uuid",
  "property_name": "Glass Transition Temperature",
  "numeric_value": -45.2,
  "text_value": null,
  "boolean_value": null,
  "confidence": 0.85,
  "calculation_method": "ML Model v1",
  "created_at": "2025-04-18T15:00:00Z"
}
```

#### Update Prediction

```
PUT /mixtures/{mixture_id}/predictions/{prediction_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture
- `prediction_id`: UUID of the prediction

Request Body:
```json
{
  "property_name": "Glass Transition Temperature",
  "value": -46.5,
  "confidence": 0.9,
  "calculation_method": "ML Model v1"
}
```

Response:
```json
{
  "id": "uuid",
  "mixture_id": "uuid",
  "property_type_id": "uuid",
  "property_name": "Glass Transition Temperature",
  "numeric_value": -46.5,
  "text_value": null,
  "boolean_value": null,
  "confidence": 0.9,
  "calculation_method": "ML Model v1",
  "created_at": "2025-04-18T15:00:00Z"
}
```

#### Delete Prediction

```
DELETE /mixtures/{mixture_id}/predictions/{prediction_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture
- `prediction_id`: UUID of the prediction

Response:
```json
{
  "message": "Prediction with ID {prediction_id} deleted successfully"
}
```

### Experiments

#### List Experiments for a Mixture

```
GET /mixtures/{mixture_id}/experiments
```

Path Parameters:
- `mixture_id`: UUID of the mixture

Response:
```json
[
  {
    "id": "uuid",
    "mixture_id": "uuid",
    "property_type_id": "uuid",
    "property_name": "Glass Transition Temperature",
    "numeric_value": -47.3,
    "text_value": null,
    "boolean_value": null,
    "experimental_conditions": "Measured using DSC at 10°C/min",
    "date_performed": "2025-04-15",
    "created_at": "2025-04-18T15:00:00Z"
  }
]
```

#### Get Experiment

```
GET /mixtures/{mixture_id}/experiments/{experiment_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture
- `experiment_id`: UUID of the experiment

Response:
```json
{
  "id": "uuid",
  "mixture_id": "uuid",
  "property_type_id": "uuid",
  "property_name": "Glass Transition Temperature",
  "numeric_value": -47.3,
  "text_value": null,
  "boolean_value": null,
  "experimental_conditions": "Measured using DSC at 10°C/min",
  "date_performed": "2025-04-15",
  "created_at": "2025-04-18T15:00:00Z"
}
```

#### Create Experiment

```
POST /mixtures/{mixture_id}/experiments
```

Path Parameters:
- `mixture_id`: UUID of the mixture

Request Body:
```json
{
  "property_name": "Glass Transition Temperature",
  "value": -47.3,
  "experimental_conditions": "Measured using DSC at 10°C/min",
  "date_performed": "2025-04-15"
}
```

Response:
```json
{
  "id": "uuid",
  "mixture_id": "uuid",
  "property_type_id": "uuid",
  "property_name": "Glass Transition Temperature",
  "numeric_value": -47.3,
  "text_value": null,
  "boolean_value": null,
  "experimental_conditions": "Measured using DSC at 10°C/min",
  "date_performed": "2025-04-15",
  "created_at": "2025-04-18T15:00:00Z"
}
```

#### Update Experiment

```
PUT /mixtures/{mixture_id}/experiments/{experiment_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture
- `experiment_id`: UUID of the experiment

Request Body:
```json
{
  "property_name": "Glass Transition Temperature",
  "value": -48.1,
  "experimental_conditions": "Measured using DSC at 5°C/min",
  "date_performed": "2025-04-16"
}
```

Response:
```json
{
  "id": "uuid",
  "mixture_id": "uuid",
  "property_type_id": "uuid",
  "property_name": "Glass Transition Temperature",
  "numeric_value": -48.1,
  "text_value": null,
  "boolean_value": null,
  "experimental_conditions": "Measured using DSC at 5°C/min",
  "date_performed": "2025-04-16",
  "created_at": "2025-04-18T15:00:00Z"
}
```

#### Delete Experiment

```
DELETE /mixtures/{mixture_id}/experiments/{experiment_id}
```

Path Parameters:
- `mixture_id`: UUID of the mixture
- `experiment_id`: UUID of the experiment

Response:
```json
{
  "message": "Experiment with ID {experiment_id} deleted successfully"
}
```

### Comparisons

#### Compare Prediction with Experiment

```
GET /mixtures/{mixture_id}/compare
```

Path Parameters:
- `mixture_id`: UUID of the mixture

Query Parameters:
- `property_name`: Name of the property to compare

Response:
```json
{
  "prediction": {
    "id": "uuid",
    "mixture_id": "uuid",
    "property_type_id": "uuid",
    "property_name": "Glass Transition Temperature",
    "numeric_value": -45.2,
    "text_value": null,
    "boolean_value": null,
    "confidence": 0.85,
    "calculation_method": "ML Model v1",
    "created_at": "2025-04-18T15:00:00Z"
  },
  "experiment": {
    "id": "uuid",
    "mixture_id": "uuid",
    "property_type_id": "uuid",
    "property_name": "Glass Transition Temperature",
    "numeric_value": -47.3,
    "text_value": null,
    "boolean_value": null,
    "experimental_conditions": "Measured using DSC at 10°C/min",
    "date_performed": "2025-04-15",
    "created_at": "2025-04-18T15:00:00Z"
  },
  "difference": 2.1,
  "percent_error": 4.44
}
```

## Authentication Endpoints

### Login

```
POST /auth/login
```

Request Body:
```json
{
  "email": "user@example.com",
  "password": "password"
}
```

Response:
```json
{
  "message": "Authentication successful",
  "user": {
    "id": "user_id",
    "email": "user@example.com"
  }
}
```

### Register

```
POST /auth/register
```

Request Body:
```json
{
  "email": "user@example.com",
  "password": "password"
}
```

Response:
```json
{
  "message": "Registration successful",
  "user": {
    "id": "user_id",
    "email": "user@example.com"
  }
}
```

### Logout

```
POST /auth/logout
```

Response:
```json
{
  "message": "Logout successful"
}
```

### Reset Password

```
POST /auth/reset-password
```

Request Body:
```json
{
  "email": "user@example.com"
}
```

Response:
```json
{
  "message": "Password reset email sent"
}
```

### Update Password

```
POST /auth/update-password
```

Request Body:
```json
{
  "password": "new_password"
}
```

Response:
```json
{
  "message": "Password updated successfully"
}
```

### Update Profile

```
POST /auth/update-profile
```

Request Body:
```json
{
  "user_data": {
    "display_name": "John Doe",
    "organization": "Research Lab"
  }
}
```

Response:
```json
{
  "message": "Profile updated successfully",
  "user": {
    "id": "user_id",
    "email": "user@example.com",
    "user_metadata": {
      "display_name": "John Doe",
      "organization": "Research Lab"
    }
  }
}