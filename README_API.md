# CryoProtect Analyzer API

This document provides information about the CryoProtect Analyzer API, which allows access to and manipulation of data in the Supabase database.

## Overview

The CryoProtect Analyzer API is a RESTful API built with Flask that provides endpoints for:

- Retrieving molecules with their properties
- Retrieving mixtures with their components
- Creating new mixtures
- Adding predictions for mixtures
- Recording experimental results
- Comparing predictions with experimental results

The API uses Supabase for authentication and data storage.

## Getting Started

### Prerequisites

- Python 3.6+
- Supabase project with the CryoProtect schema applied
- Required Python packages (see `requirements.txt`)

### Installation

1. Clone the repository
2. Install the required packages:
   ```
   pip install -r requirements.txt
   ```
3. Create a `.env` file with your Supabase credentials:
   ```
   SUPABASE_URL=https://your-project-ref.supabase.co
   SUPABASE_KEY=your-anon-key
   SUPABASE_USER=user@example.com
   SUPABASE_PASSWORD=your-password
   ```
4. Run the API:
   ```
   python app.py
   ```

## API Documentation

The API documentation is available at `/swagger-ui/` when the API is running. This provides a Swagger UI interface for exploring and testing the API endpoints.

### Authentication

The API uses Supabase authentication. To authenticate, send a POST request to `/auth/login`. The API will use the credentials from the `.env` file.

For protected endpoints, include the JWT token in the Authorization header:
```
Authorization: Bearer <token>
```

### Endpoints

#### Molecules

- `GET /api/v1/molecules` - Get a list of all molecules with their properties
- `POST /api/v1/molecules` - Import a molecule from PubChem (requires authentication)
- `GET /api/v1/molecules/{molecule_id}` - Get a specific molecule with its properties

#### Mixtures

- `GET /api/v1/mixtures` - Get a list of all mixtures with their components
- `POST /api/v1/mixtures` - Create a new mixture (requires authentication)
- `GET /api/v1/mixtures/{mixture_id}` - Get a specific mixture with its components

#### Predictions

- `GET /api/v1/mixtures/{mixture_id}/predictions` - Get all predictions for a mixture
- `POST /api/v1/mixtures/{mixture_id}/predictions` - Add a prediction for a mixture (requires authentication)

#### Experiments

- `GET /api/v1/mixtures/{mixture_id}/experiments` - Get all experiments for a mixture
- `POST /api/v1/mixtures/{mixture_id}/experiments` - Record an experiment for a mixture (requires authentication)

#### Comparisons

- `GET /api/v1/mixtures/{mixture_id}/comparisons?property_name={property_name}` - Compare prediction with experiment for a mixture and property

## Request and Response Examples

### Molecules

#### Get all molecules

```
GET /api/v1/molecules
```

Response:
```json
[
  {
    "id": "123e4567-e89b-12d3-a456-426614174000",
    "cid": 753,
    "name": "Glycerol",
    "molecular_formula": "C3H8O3",
    "smiles": "C(C(CO)O)O",
    "pubchem_link": "https://pubchem.ncbi.nlm.nih.gov/compound/753",
    "created_at": "2025-04-14T21:00:00Z",
    "updated_at": "2025-04-14T21:00:00Z",
    "properties": {
      "Molecular Weight": 92.09,
      "LogP": -1.76,
      "TPSA": 60.69,
      "H-Bond Donors": 3,
      "H-Bond Acceptors": 3,
      "Toxicity": "Low toxicity",
      "Stability": "Stable under normal conditions",
      "Environmental Safety": "Biodegradable",
      "Total Score": 180
    }
  }
]
```

#### Import a molecule from PubChem

```
POST /api/v1/molecules
Content-Type: application/json
Authorization: Bearer <token>

{
  "cid": 753
}
```

Response:
```json
{
  "id": "123e4567-e89b-12d3-a456-426614174000",
  "cid": 753,
  "name": "Glycerol",
  "molecular_formula": "C3H8O3",
  "smiles": "C(C(CO)O)O",
  "pubchem_link": "https://pubchem.ncbi.nlm.nih.gov/compound/753",
  "created_at": "2025-04-14T21:00:00Z",
  "updated_at": "2025-04-14T21:00:00Z",
  "properties": {
    "Molecular Weight": 92.09,
    "LogP": -1.76,
    "TPSA": 60.69,
    "H-Bond Donors": 3,
    "H-Bond Acceptors": 3,
    "Toxicity": "Low toxicity",
    "Stability": "Stable under normal conditions",
    "Environmental Safety": "Biodegradable",
    "Total Score": 180
  }
}
```

### Mixtures

#### Create a new mixture

```
POST /api/v1/mixtures
Content-Type: application/json
Authorization: Bearer <token>

{
  "name": "Glycerol-Water Solution",
  "description": "A 30% glycerol solution in water",
  "components": [
    {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
      "concentration": 30,
      "concentration_unit": "%"
    },
    {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174001",
      "concentration": 70,
      "concentration_unit": "%"
    }
  ]
}
```

Response:
```json
{
  "id": "123e4567-e89b-12d3-a456-426614174002",
  "name": "Glycerol-Water Solution",
  "description": "A 30% glycerol solution in water",
  "created_at": "2025-04-14T21:00:00Z",
  "updated_at": "2025-04-14T21:00:00Z",
  "components": [
    {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
      "cid": 753,
      "name": "Glycerol",
      "concentration": 30,
      "concentration_unit": "%"
    },
    {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174001",
      "cid": 962,
      "name": "Water",
      "concentration": 70,
      "concentration_unit": "%"
    }
  ]
}
```

### Predictions

#### Add a prediction

```
POST /api/v1/mixtures/123e4567-e89b-12d3-a456-426614174002/predictions
Content-Type: application/json
Authorization: Bearer <token>

{
  "property_name": "Freezing Point",
  "value": -15.3,
  "confidence": 0.9,
  "calculation_method": "CryoProtect Scoring"
}
```

Response:
```json
{
  "id": "123e4567-e89b-12d3-a456-426614174003",
  "mixture_id": "123e4567-e89b-12d3-a456-426614174002",
  "property_type_id": "123e4567-e89b-12d3-a456-426614174004",
  "property_name": "Freezing Point",
  "numeric_value": -15.3,
  "text_value": null,
  "boolean_value": null,
  "confidence": 0.9,
  "calculation_method": "CryoProtect Scoring",
  "created_at": "2025-04-14T21:00:00Z"
}
```

### Experiments

#### Record an experiment

```
POST /api/v1/mixtures/123e4567-e89b-12d3-a456-426614174002/experiments
Content-Type: application/json
Authorization: Bearer <token>

{
  "property_name": "Freezing Point",
  "value": -14.8,
  "experimental_conditions": "Standard pressure, cooling rate 1°C/min",
  "date_performed": "2025-04-14"
}
```

Response:
```json
{
  "id": "123e4567-e89b-12d3-a456-426614174005",
  "mixture_id": "123e4567-e89b-12d3-a456-426614174002",
  "property_type_id": "123e4567-e89b-12d3-a456-426614174004",
  "property_name": "Freezing Point",
  "numeric_value": -14.8,
  "text_value": null,
  "boolean_value": null,
  "experimental_conditions": "Standard pressure, cooling rate 1°C/min",
  "date_performed": "2025-04-14",
  "created_at": "2025-04-14T21:00:00Z"
}
```

### Comparisons

#### Compare prediction with experiment

```
GET /api/v1/mixtures/123e4567-e89b-12d3-a456-426614174002/comparisons?property_name=Freezing%20Point
```

Response:
```json
{
  "prediction": {
    "numeric_value": -15.3,
    "text_value": null,
    "boolean_value": null,
    "confidence": 0.9,
    "method": "CryoProtect Scoring"
  },
  "experiment": {
    "numeric_value": -14.8,
    "text_value": null,
    "boolean_value": null,
    "conditions": "Standard pressure, cooling rate 1°C/min",
    "date": "2025-04-14"
  },
  "difference": 0.5,
  "percent_error": 3.38
}
```

## Error Handling

The API returns appropriate HTTP status codes and error messages for different error conditions:

- 400 Bad Request - Invalid request data
- 401 Unauthorized - Missing or invalid authentication token
- 403 Forbidden - Insufficient permissions
- 404 Not Found - Resource not found
- 409 Conflict - Resource already exists
- 500 Internal Server Error - Server-side error

Example error response:
```json
{
  "message": "Resource not found"
}
```

## Development

### Running Tests

```
python -m unittest discover tests
```

### Adding New Endpoints

1. Define the resource class in `api/resources.py`
2. Register the resource in `api/__init__.py`
3. Update the API documentation

## License

This project is licensed under the MIT License - see the LICENSE file for details.