# CryoProtect Analyzer API

...

## Toxicity Data API

### `GET /api/toxicity/molecule/{molecule_id}`

Retrieves toxicity data for a specific molecule from the Tox21 database.

**Parameters:**
- `molecule_id` (path parameter): UUID of the molecule
- `assay_id` (query parameter, optional): Filter by assay ID
- `endpoint` (query parameter, optional): Filter by toxicological endpoint
- `active_only` (query parameter, optional): Return only active results (default: false)

**Example Request:**
```http
GET /api/toxicity/molecule/123e4567-e89b-12d3-a456-426614174000
```

**Example Response:**
```json
{
  "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
  "toxicity_data": [
    {
      "id": "789e0123-c45d-67e8-f901-234567890123",
      "assay_id": "456e7890-a12b-34c5-d67e-890123456789",
      "activity_value": 10.5,
      "activity_unit": "µM",
      "hit_call": true,
      "reliability_score": 0.85
    }
  ],
  "count": 1
}
```

### `GET /api/toxicity/endpoints`

Retrieves available toxicity endpoints from the Tox21 database.

**Example Request:**
```http
GET /api/toxicity/endpoints
```

**Example Response:**
```json
{
  "endpoints": [
    {
      "id": "123e4567-e89b-12d3-a456-426614174000",
      "name": "nuclear_receptor",
      "description": "Nuclear receptor activity",
      "category": "endocrine"
    },
    {
      "id": "234e5678-f90a-12b3-c456-789012345678",
      "name": "stress_response",
      "description": "Cellular stress response pathway",
      "category": "cellular"
    }
  ],
  "count": 2
}
```

### `GET /api/toxicity/assays`

Retrieves available toxicity assays from the Tox21 database.

**Parameters:**
- `endpoint` (query parameter, optional): Filter by toxicological endpoint

**Example Request:**
```http
GET /api/toxicity/assays?endpoint=nuclear_receptor
```

**Example Response:**
```json
{
  "assays": [
    {
      "id": "456e7890-a12b-34c5-d67e-890123456789",
      "name": "NR-AR",
      "description": "Androgen Receptor",
      "toxicological_endpoint": "nuclear_receptor",
      "assay_source": "Tox21"
    },
    {
      "id": "567e8901-b23c-45d6-e78f-901234567890",
      "name": "NR-ER",
      "description": "Estrogen Receptor",
      "toxicological_endpoint": "nuclear_receptor",
      "assay_source": "Tox21"
    }
  ],
  "count": 2
}
```

### `GET /api/toxicity/scores/molecule/{molecule_id}`

Retrieves toxicity scores for a specific molecule. If scores don't exist, they will be calculated automatically.

**Parameters:**
- `molecule_id` (path parameter): UUID of the molecule

**Example Request:**
```http
GET /api/toxicity/scores/molecule/123e4567-e89b-12d3-a456-426614174000
```

**Example Response:**
```json
{
  "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
  "scores": [
    {
      "id": "890e1234-d56e-78f9-0123-456789012345",
      "endpoint_id": "123e4567-e89b-12d3-a456-426614174000",
      "score_type": "toxicity_score",
      "score_value": 0.35,
      "confidence": 0.85,
      "endpoint": {
        "id": "123e4567-e89b-12d3-a456-426614174000",
        "name": "nuclear_receptor",
        "description": "Nuclear receptor activity",
        "category": "endocrine"
      }
    }
  ],
  "count": 1
}
```

### `GET /api/toxicity/scores/mixture/{mixture_id}`

Retrieves toxicity scores for a specific mixture, including both aggregate scores and component scores.

**Parameters:**
- `mixture_id` (path parameter): UUID of the mixture

**Example Request:**
```http
GET /api/toxicity/scores/mixture/345e6789-0a1b-23c4-d56e-789012345678
```

**Example Response:**
```json
{
  "mixture_id": "345e6789-0a1b-23c4-d56e-789012345678",
  "aggregate_score": {
    "overall_toxicity": 0.42,
    "confidence": 0.75
  },
  "component_scores": [
    {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
      "name": "DMSO",
      "concentration": 0.6,
      "toxicity_score": 0.35
    },
    {
      "molecule_id": "234e5678-f90a-12b3-c456-789012345678",
      "name": "Glycerol",
      "concentration": 0.4,
      "toxicity_score": 0.25
    }
  ],
  "component_count": 2
}
```

### `GET /api/toxicity/unified/molecule/{molecule_id}`

Retrieves a unified score for a molecule that combines efficacy, toxicity, and glass transition temperature (Tg) data.

**Parameters:**
- `molecule_id` (path parameter): UUID of the molecule
- `context` (query parameter, optional): Application context for weighting (default: "general")
- `algorithm` (query parameter, optional): Algorithm to use for predictive models (default: "random_forest")
- `recalculate` (query parameter, optional): Whether to recalculate scores even if they already exist (default: false)

**Example Request:**
```http
GET /api/toxicity/unified/molecule/123e4567-e89b-12d3-a456-426614174000?context=cell_preservation
```

**Example Response:**
```json
{
  "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
  "name": "DMSO",
  "unified_score": 78.5,
  "application_context": "cell_preservation",
  "component_scores": {
    "efficacy": {
      "score": 85.0,
      "weight": 0.5,
      "weighted_score": 42.5
    },
    "toxicity": {
      "score": 65.0,
      "weight": 0.25,
      "weighted_score": 16.25
    },
    "glass_transition": {
      "score": 79.0,
      "weight": 0.25,
      "weighted_score": 19.75
    }
  },
  "timestamp": "2025-04-22T23:45:12Z"
}
```

### `GET /api/toxicity/unified/mixture/{mixture_id}`

Retrieves a unified score for a mixture that combines efficacy, toxicity, and glass transition temperature (Tg) data.

**Parameters:**
- `mixture_id` (path parameter): UUID of the mixture
- `context` (query parameter, optional): Application context for weighting (default: "general")
- `algorithm` (query parameter, optional): Algorithm to use for predictive models (default: "random_forest")
- `recalculate` (query parameter, optional): Whether to recalculate scores even if they already exist (default: false)

**Example Request:**
```http
GET /api/toxicity/unified/mixture/345e6789-0a1b-23c4-d56e-789012345678?context=organ_preservation
```

**Example Response:**
```json
{
  "mixture_id": "345e6789-0a1b-23c4-d56e-789012345678",
  "name": "DMSO-Glycerol Mixture",
  "unified_score": 82.3,
  "application_context": "organ_preservation",
  "component_scores": {
    "efficacy": {
      "score": 80.0,
      "weight": 0.3,
      "weighted_score": 24.0
    },
    "toxicity": {
      "score": 88.0,
      "weight": 0.4,
      "weighted_score": 35.2
    },
    "glass_transition": {
      "score": 77.0,
      "weight": 0.3,
      "weighted_score": 23.1
    }
  },
  "components": [
    {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
      "name": "DMSO",
      "concentration": 0.6
    },
    {
      "molecule_id": "234e5678-f90a-12b3-c456-789012345678",
      "name": "Glycerol",
      "concentration": 0.4
    }
  ],
  "timestamp": "2025-04-22T23:46:05Z"
}
```

### `GET /api/toxicity/unified/contexts`

Retrieves available application contexts for unified scoring, including their descriptions and weight profiles.

**Example Request:**
```http
GET /api/toxicity/unified/contexts
```

**Example Response:**
```json
{
  "contexts": {
    "general": {
      "description": "Balanced profile suitable for general cryopreservation applications",
      "weights": {
        "efficacy": 0.4,
        "toxicity": 0.3,
        "glass_transition": 0.3
      },
      "tg_optimal_range": [-120, -80]
    },
    "cell_preservation": {
      "description": "Prioritizes efficacy and permeability for cell preservation",
      "weights": {
        "efficacy": 0.5,
        "toxicity": 0.25,
        "glass_transition": 0.25
      },
      "tg_optimal_range": [-130, -100]
    },
    "organ_preservation": {
      "description": "Prioritizes low toxicity and optimal glass transition for organ preservation",
      "weights": {
        "efficacy": 0.3,
        "toxicity": 0.4,
        "glass_transition": 0.3
      },
      "tg_optimal_range": [-110, -90]
    },
    "long_term_storage": {
      "description": "Prioritizes glass-forming ability and stability for long-term storage",
      "weights": {
        "efficacy": 0.3,
        "toxicity": 0.2,
        "glass_transition": 0.5
      },
      "tg_optimal_range": [-100, -70]
    },
    "sensitive_tissues": {
      "description": "Prioritizes safety for sensitive tissues like neural or reproductive cells",
      "weights": {
        "efficacy": 0.3,
        "toxicity": 0.5,
        "glass_transition": 0.2
      },
      "tg_optimal_range": [-120, -100]
    }
  }
}
```

### `POST /api/toxicity/unified/batch`

Calculates unified scores for multiple molecules or mixtures in a single request.

**Request Body:**
```json
{
  "entity_ids": ["id1", "id2", "id3"],
  "entity_type": "molecule", // or "mixture"
  "application_context": "general", // optional
  "algorithm": "random_forest", // optional
  "recalculate": false // optional
}
```

**Example Request:**
```http
POST /api/toxicity/unified/batch
Content-Type: application/json

{
  "entity_ids": [
    "123e4567-e89b-12d3-a456-426614174000",
    "234e5678-f90a-12b3-c456-789012345678"
  ],
  "entity_type": "molecule",
  "application_context": "cell_preservation"
}
```

**Example Response:**
```json
{
  "results": {
    "123e4567-e89b-12d3-a456-426614174000": {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
      "name": "DMSO",
      "unified_score": 78.5,
      "application_context": "cell_preservation",
      "component_scores": {
        "efficacy": {
          "score": 85.0,
          "weight": 0.5,
          "weighted_score": 42.5
        },
        "toxicity": {
          "score": 65.0,
          "weight": 0.25,
          "weighted_score": 16.25
        },
        "glass_transition": {
          "score": 79.0,
          "weight": 0.25,
          "weighted_score": 19.75
        }
      },
      "timestamp": "2025-04-22T23:47:30Z"
    },
    "234e5678-f90a-12b3-c456-789012345678": {
      "molecule_id": "234e5678-f90a-12b3-c456-789012345678",
      "name": "Glycerol",
      "unified_score": 72.8,
      "application_context": "cell_preservation",
      "component_scores": {
        "efficacy": {
          "score": 75.0,
          "weight": 0.5,
          "weighted_score": 37.5
        },
        "toxicity": {
          "score": 80.0,
          "weight": 0.25,
          "weighted_score": 20.0
        },
        "glass_transition": {
          "score": 61.2,
          "weight": 0.25,
          "weighted_score": 15.3
        }
      },
      "timestamp": "2025-04-22T23:47:30Z"
    }
  },
  "errors": [],
  "summary": {
    "total": 2,
    "successful": 2,
    "failed": 0,
    "application_context": "cell_preservation"
  }
}
```

## Compare Molecular and Mixture Properties

### `POST /api/v1/compare-properties`

Compare key properties of selected molecules and/or mixtures side-by-side. This endpoint is useful for visualizing differences in molecular weight, logP, TPSA, toxicity, vitrification probability, and other scientific properties.

**Request Body:**
```json
{
  "ids": [
    "molecule_id_1",
    "mixture_id_2",
    "molecule_id_3"
  ]
}
```
- `ids`: List of molecule and/or mixture IDs to compare. IDs can be mixed.

**Example Request:**
```http
POST /api/v1/compare-properties
Content-Type: application/json

{
  "ids": ["mol_abc123", "mix_xyz789", "mol_def456"]
}
```

**Example Response:**
```json
{
  "comparison": [
    {
      "id": "mol_abc123",
      "name": "Ethylene Glycol",
      "type": "molecule",
      "molecular_weight": 62.07,
      "logP": -1.36,
      "TPSA": 40.46,
      "toxicity": "Low",
      "vitrification_probability": 0.85
    },
    {
      "id": "mix_xyz789",
      "name": "EG + DMSO Mixture",
      "type": "mixture",
      "molecular_weight": 90.12,
      "logP": -0.95,
      "TPSA": 60.00,
      "toxicity": "Moderate",
      "vitrification_probability": 0.92
    },
    {
      "id": "mol_def456",
      "name": "DMSO",
      "type": "molecule",
      ...
    }
  ]
}
```

---

## Batch Operations

### `POST /api/v1/batch`

Perform batch operations on molecules, mixtures, or experiments. Supported operations include property calculation, mixture optimization, predictive scoring, and export of protocols/results.

**Request Body:**
```json
{
  "operation": "property_calculation", // or "mixture_optimization", "predictive_scoring", "export"
  "entity_type": "molecule", // or "mixture", "experiment"
  "ids": ["id1", "id2", "id3"]
}
```
- `operation`: The batch operation to perform.
- `entity_type`: The type of entity to process.
- `ids`: List of IDs to process in batch.

**Example Request:**
```http
POST /api/v1/batch
Content-Type: application/json

{
  "operation": "property_calculation",
  "entity_type": "mixture",
  "ids": ["mix_001", "mix_002", "mix_003"]
}
```

**Example Response:**
```json
{
  "status": "SUCCESS",
  "results": [
    {"id": "mix_001", "result": { /* property data */ }},
    {"id": "mix_002", "result": { /* property data */ }},
    {"id": "mix_003", "result": { /* property data */ }}
  ],
  "errors": []
}
```

**Partial Failure Example:**
```json
{
  "status": "COMPLETED_WITH_WARNINGS",
  "results": [
    {"id": "mix_001", "result": { /* property data */ }},
    {"id": "mix_003", "result": { /* property data */ }}
  ],
  "errors": [
    {"id": "mix_002", "error": "Mixture not found"}
  ]
}
```

**Error Example:**
```json
{
  "status": "ERROR",
  "results": [],
  "errors": [
    {"id": "mix_002", "error": "Mixture not found"},
    {"id": "mix_004", "error": "Invalid ID"}
  ]
}
```

**Supported Operations:**
- `property_calculation`: Calculate properties for molecules or mixtures.
- `mixture_optimization`: Optimize mixtures (entity_type must be "mixture").
- `predictive_scoring`: Predict scores for molecules, mixtures, or experiments.
- `export`: Export protocols/results for the given entities.

**Notes:**
- The endpoint returns a list of results and a list of errors for partial failures.
- For export operations, the result may include file links or protocol data.
- All operations require authentication (token).
- For frontend integration, use the `status` field to determine if the operation was fully or partially successful.

---

## Data Export: Protocols and Results

### `POST /api/v1/export`

Export protocols, experiment results, molecules, mixtures, predictions, or comparisons in CSV, JSON, Excel, or PDF format.

**Request Body:**
```json
{
  "format": "csv",         // "csv", "json", "excel", or "pdf"
  "data_type": "protocols", // "protocols", "experiments", "molecules", "mixtures", "predictions", or "comparisons"
  "id": "protocol_id_123", // Optional: ID of the item to export (required for protocols)
  "include_related": false // Optional: Include related data (where supported)
}
```
- `format`: Export file format.
- `data_type`: Type of data to export.
- `id`: ID of the protocol, experiment, etc. (required for protocols).
- `include_related`: Include related predictions/experiments (for mixtures).

**Example: Export a Protocol as JSON**
```http
POST /api/v1/export
Content-Type: application/json

{
  "format": "json",
  "data_type": "protocols",
  "id": "protocol_id_123"
}
```

**Example: Export Experiment Results as CSV**
```http
POST /api/v1/export
Content-Type: application/json

{
  "format": "csv",
  "data_type": "experiments",
  "id": "mixture_id_456"
}
```

**Response:**
- Returns a downloadable file (CSV, JSON, Excel, or PDF) or an error message.

**Supported data_type values:**
- `molecules`
- `mixtures`
- `predictions`
- `experiments` (experiment results)
- `comparisons`
- `protocols` (stepwise protocols, by protocol ID)

**Notes:**
- Protocol export requires a valid protocol ID. Protocol export depends on the implementation of `ProtocolDesigner.get_saved_protocol`. If protocol storage is not implemented, export may return an error.
- For experiment results, use `data_type: "experiments"` and provide the mixture ID.
- All export operations require authentication (token).
- For frontend integration, use the returned file or error message to inform the user.
