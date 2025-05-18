# Experimental Data API Reference

This document provides comprehensive documentation for the Enhanced Experimental Data System API. It covers all available endpoints, their parameters, response formats, and usage examples.

## Table of Contents

- [Convex API](#convex-api)
  - [Authentication](#authentication)
  - [Response Format](#response-format)
  - [Experiments](#experiments)
  - [Experiment Results](#experiment-results)
  - [Protocols](#protocols)
  - [Tissue Types](#tissue-types)
  - [Equipment](#equipment)
  - [Time Series](#time-series)
  - [Validation Rules](#validation-rules)
- [HTTP API](#http-api)
  - [Authentication](#http-authentication)
  - [Response Format](#http-response-format)
  - [Endpoints](#http-endpoints)
- [Error Handling](#error-handling)
- [Pagination and Filtering](#pagination-and-filtering)

## Convex API

The primary API for the Enhanced Experimental Data System is implemented using Convex. This API provides direct access to the database and business logic for frontend applications.

### Authentication

The API uses Convex's built-in authentication system. All operations require an authenticated user and respect access control rules:

- Public experiments are accessible to all users
- Private experiments are only accessible to the creator
- Some operations (like delete) are restricted to the experiment creator

### Response Format

Convex API functions return structured TypeScript objects with proper type definitions. All responses include:

- Core data fields
- References to related entities
- Timestamp information
- Access control metadata

Errors are thrown as standard JavaScript exceptions with descriptive messages.

### Experiments

#### Create Experiment

Creates a new enhanced experiment record.

```typescript
// Function: experiments.createEnhancedExperiment
// Input:
const experimentId = await client.mutation("experiments.createEnhancedExperiment", {
  experiment: {
    name: "DMSO Vitrification Test",
    description: "Testing vitrification with DMSO concentration series",
    experimentTypeId: "vitrification",
    protocolId: "conv:protocols:abc123", // Protocol ID
    mixtureId: "conv:mixtures:def456",   // Mixture ID
    temperature: -196,
    temperatureUnit: "C",
    coolingRate: -200,
    coolingRateUnit: "C/min",
    thawingRate: 10,
    thawingRateUnit: "C/min",
    pressure: 101.3,
    pressureUnit: "kPa",
    parameters: {
      sampleCount: 10,
      controlGroup: true
    },
    provenance: {
      method: "Standard Operating Procedure v2.3",
      operator: "Dr. Smith"
    },
    projectId: "conv:projects:ghi789", // Optional project ID
    date: Date.now(),
    status: "planned", // planned, in-progress, completed, failed
    tags: ["vitrification", "dmso", "mouse-embryos"],
    public: false // Default is false
  }
});
```

#### Get Experiment

Retrieves an experiment by ID with optional related data.

```typescript
// Function: experiments.getEnhancedExperiment
// Input:
const experiment = await client.query("experiments.getEnhancedExperiment", {
  experimentId: "conv:enhancedExperiments:abc123",
  options: {
    includeResults: true,
    includeProtocol: true,
    includeMixture: true,
    includeTissueTypes: true,
    includeEquipment: true,
    includeTimeSeries: true
  }
});
```

#### Update Experiment

Updates an existing experiment.

```typescript
// Function: experiments.updateEnhancedExperiment
// Input:
await client.mutation("experiments.updateEnhancedExperiment", {
  experimentId: "conv:enhancedExperiments:abc123",
  update: {
    name: "Updated DMSO Test",
    description: "Updated description",
    status: "in-progress",
    tags: ["updated", "dmso", "mouse-embryos"]
  }
});
```

#### Delete Experiment

Deletes an experiment and all its related data (results, equipment links, time series).

```typescript
// Function: experiments.deleteEnhancedExperiment
// Input:
await client.mutation("experiments.deleteEnhancedExperiment", {
  experimentId: "conv:enhancedExperiments:abc123"
});
```

#### List Experiments

Retrieves a list of experiments with filtering, sorting, and pagination.

```typescript
// Function: experiments.listEnhancedExperiments
// Input:
const experiments = await client.query("experiments.listEnhancedExperiments", {
  filter: {
    name: "DMSO", // Optional name filter (partial match)
    experimentTypeId: "vitrification", // Optional type filter
    protocolId: "conv:protocols:abc123", // Optional protocol filter
    mixtureId: "conv:mixtures:def456", // Optional mixture filter
    conductedBy: "conv:users:ghi789", // Optional user filter
    projectId: "conv:projects:jkl012", // Optional project filter
    status: "completed", // Optional status filter
    dateRange: {
      start: new Date('2025-01-01').getTime(),
      end: new Date('2025-12-31').getTime()
    },
    tags: ["dmso", "vitrification"], // Optional tag filter (experiments with any of these tags)
    public: true, // Optional visibility filter
    tissueTypeId: "conv:tissueTypes:mno345" // Optional tissue type filter
  },
  options: {
    limit: 20, // Number of experiments to return
    cursor: "cursor_token", // Pagination cursor from previous query
    includeResults: true, // Include experiment results
    includeProtocol: true, // Include protocol details
    includeMixture: true, // Include mixture details
    includeTissueTypes: true, // Include tissue type details
    includeEquipment: true, // Include equipment details
    includeTimeSeries: true, // Include time series details
    sortBy: "date", // Sort by field (name, date, status, updatedAt)
    sortDirection: "desc" // Sort direction (asc, desc)
  }
});
```

#### Search Experiments

Searches for experiments by name.

```typescript
// Function: experiments.searchEnhancedExperiments
// Input:
const experiments = await client.query("experiments.searchEnhancedExperiments", {
  query: "DMSO vitrification",
  limit: 20
});
```

#### Update Experiment Status

Updates only the status of an experiment.

```typescript
// Function: experiments.updateEnhancedExperimentStatus
// Input:
await client.mutation("experiments.updateEnhancedExperimentStatus", {
  experimentId: "conv:enhancedExperiments:abc123",
  status: "completed" // planned, in-progress, completed, failed
});
```

### Experiment Results

#### Create Result

Adds a new result to an experiment.

```typescript
// Function: experiments.createEnhancedExperimentResult
// Input:
const resultId = await client.mutation("experiments.createEnhancedExperimentResult", {
  result: {
    experimentId: "conv:enhancedExperiments:abc123",
    moleculeId: "conv:molecules:def456", // Optional molecule reference
    mixtureId: "conv:mixtures:ghi789", // Optional mixture reference
    tissueTypeId: "conv:tissueTypes:jkl012", // Tissue type reference
    parameterName: "viability_percentage",
    value: 87.5, // Can be string, number, boolean, or null
    numericValue: 87.5, // Numeric representation for filtering/sorting
    units: "%",
    uncertainty: {
      type: "standard_deviation", // standard_deviation, range, confidence_interval
      value: 4.2, // Single value for std dev, array [min, max] for range/CI
      confidence: 0.95 // For confidence intervals (e.g., 0.95 for 95% CI)
    },
    provenance: {
      method: "Flow Cytometry",
      reference: "Standard Operating Procedure v2.3",
      timestamp: Date.now(),
      operator: "Dr. Smith"
    },
    notes: "Analyzed using fluorescent markers"
  }
});
```

#### Get Results for Experiment

Retrieves all results for a specific experiment.

```typescript
// Function: experiments.getEnhancedExperimentResults
// Input:
const results = await client.query("experiments.getEnhancedExperimentResults", {
  experimentId: "conv:enhancedExperiments:abc123"
});
```

#### Update Result

Updates an existing result.

```typescript
// Function: experiments.updateEnhancedExperimentResult
// Input:
await client.mutation("experiments.updateEnhancedExperimentResult", {
  resultId: "conv:enhancedExperimentResults:abc123",
  update: {
    value: 88.2,
    numericValue: 88.2,
    uncertainty: {
      type: "standard_deviation",
      value: 3.8
    },
    notes: "Updated after recalibration"
  }
});
```

#### Delete Result

Deletes a specific result.

```typescript
// Function: experiments.deleteEnhancedExperimentResult
// Input:
await client.mutation("experiments.deleteEnhancedExperimentResult", {
  resultId: "conv:enhancedExperimentResults:abc123"
});
```

### Protocols

#### Create Protocol

Creates a new protocol template or instance.

```typescript
// Function: experiments.createProtocol
// Input:
const protocolId = await client.mutation("experiments.createProtocol", {
  protocol: {
    name: "Standard Vitrification Protocol",
    description: "A standard protocol for vitrification of mouse embryos",
    steps: [
      {
        name: "Equilibration Step",
        description: "Equilibrate embryos in vitrification solution",
        duration: 10,
        durationUnit: "minutes",
        temperature: 4,
        temperatureUnit: "C",
        parameters: {
          bufferSolution: "PBS",
          stirringRate: 5
        }
      },
      {
        name: "Vitrification Step",
        description: "Rapidly cool embryos to achieve vitrification",
        temperature: -196,
        temperatureUnit: "C"
      }
    ],
    parameters: {
      targetConcentration: 10.0,
      sampleType: "embryo"
    },
    isTemplate: true,
    category: "vitrification",
    public: true
  }
});
```

#### Get Protocol

Retrieves a protocol by ID.

```typescript
// Function: experiments.getProtocol
// Input:
const protocol = await client.query("experiments.getProtocol", {
  protocolId: "conv:protocols:abc123"
});
```

#### Update Protocol

Updates an existing protocol.

```typescript
// Function: experiments.updateProtocol
// Input:
await client.mutation("experiments.updateProtocol", {
  protocolId: "conv:protocols:abc123",
  update: {
    name: "Updated Vitrification Protocol",
    description: "Enhanced protocol with optimized parameters",
    steps: [
      // Updated steps
    ]
  }
});
```

#### List Protocols

Retrieves a list of protocols with filtering.

```typescript
// Function: experiments.listProtocols
// Input:
const protocols = await client.query("experiments.listProtocols", {
  filter: {
    name: "Vitrification",
    category: "cryopreservation",
    isTemplate: true,
    public: true
  },
  options: {
    limit: 20,
    cursor: "cursor_token",
    sortBy: "name",
    sortDirection: "asc"
  }
});
```

#### Create Protocol Version

Creates a new version of an existing protocol.

```typescript
// Function: experiments.createProtocolVersion
// Input:
const newVersionId = await client.mutation("experiments.createProtocolVersion", {
  parentProtocolId: "conv:protocols:abc123",
  changes: {
    name: "Vitrification Protocol v2",
    description: "Updated version with improved steps",
    steps: [
      // Modified steps
    ]
  }
});
```

### Tissue Types

#### Create Tissue Type

Creates a new tissue type record.

```typescript
// Function: experiments.createTissueType
// Input:
const tissueTypeId = await client.mutation("experiments.createTissueType", {
  tissueType: {
    name: "Mouse Embryo",
    description: "Preimplantation mouse embryo",
    species: "Mus musculus",
    taxonomyId: 10090,
    properties: {
      stage: "blastocyst",
      age: "3.5 days"
    },
    category: "embryo",
    source: "Laboratory bred",
    public: true
  }
});
```

#### Get Tissue Type

Retrieves a tissue type by ID.

```typescript
// Function: experiments.getTissueType
// Input:
const tissueType = await client.query("experiments.getTissueType", {
  tissueTypeId: "conv:tissueTypes:abc123"
});
```

#### List Tissue Types

Retrieves a list of tissue types with filtering.

```typescript
// Function: experiments.listTissueTypes
// Input:
const tissueTypes = await client.query("experiments.listTissueTypes", {
  filter: {
    name: "embryo",
    species: "Mus musculus",
    category: "embryo",
    public: true
  },
  options: {
    limit: 20,
    cursor: "cursor_token"
  }
});
```

### Equipment

#### Create Equipment Type

Creates a new equipment type category.

```typescript
// Function: experiments.createEquipmentType
// Input:
const equipmentTypeId = await client.mutation("experiments.createEquipmentType", {
  equipmentType: {
    name: "Controlled Rate Freezer",
    description: "Programmable freezing device for controlled cooling rates",
    manufacturer: "CryoMed",
    category: "freezing-equipment"
  }
});
```

#### Create Equipment

Creates a specific equipment instance.

```typescript
// Function: experiments.createEquipment
// Input:
const equipmentId = await client.mutation("experiments.createEquipment", {
  equipment: {
    equipmentTypeId: "conv:equipmentTypes:abc123",
    name: "CryoMed 7456",
    model: "7456-Advanced",
    serialNumber: "CM74560123",
    description: "Lab 3 controlled rate freezer",
    calibrationDate: new Date('2025-01-15').getTime(),
    nextCalibrationDate: new Date('2025-07-15').getTime(),
    location: "Lab 3, Bench 2",
    parameters: {
      minTemp: -180,
      maxTemp: 50,
      maxCoolingRate: 50
    }
  }
});
```

#### Link Equipment to Experiment

Associates a piece of equipment with an experiment.

```typescript
// Function: experiments.createExperimentEquipment
// Input:
const linkId = await client.mutation("experiments.createExperimentEquipment", {
  experimentEquipment: {
    experimentId: "conv:enhancedExperiments:abc123",
    equipmentId: "conv:equipment:def456",
    role: "cooling_device",
    parameters: {
      programName: "Standard Slow Cooling",
      startTemp: 20,
      endTemp: -80
    },
    notes: "Used in conjunction with LN2 plunge freezing"
  }
});
```

#### Get Equipment for Experiment

Retrieves all equipment used in a specific experiment.

```typescript
// Function: experiments.getExperimentEquipment
// Input:
const equipment = await client.query("experiments.getExperimentEquipment", {
  experimentId: "conv:enhancedExperiments:abc123"
});
```

### Time Series

#### Create Time Series

Creates a new time series definition for an experiment.

```typescript
// Function: experiments.createTimeSeries
// Input:
const timeSeriesId = await client.mutation("experiments.createTimeSeries", {
  timeSeries: {
    experimentId: "conv:enhancedExperiments:abc123",
    name: "Cooling Curve",
    description: "Temperature measurements during cooling process",
    parameterName: "temperature",
    units: "C",
    metadata: {
      sensorType: "thermocouple",
      sensorLocation: "sample-center"
    }
  }
});
```

#### Add Time Series Data Point

Adds a data point to a time series.

```typescript
// Function: experiments.addTimeSeriesDataPoint
// Input:
const dataPointId = await client.mutation("experiments.addTimeSeriesDataPoint", {
  dataPoint: {
    timeSeriesId: "conv:timeSeries:abc123",
    timestamp: Date.now(),
    value: -5.2,
    uncertainty: 0.1,
    metadata: {
      sensorReading: "stable"
    }
  }
});
```

#### Get Time Series with Data

Retrieves a time series with all its data points.

```typescript
// Function: experiments.getTimeSeriesWithData
// Input:
const timeSeriesWithData = await client.query("experiments.getTimeSeriesWithData", {
  timeSeriesId: "conv:timeSeries:abc123",
  options: {
    timeRange: {
      start: new Date('2025-05-01T10:00:00').getTime(),
      end: new Date('2025-05-01T11:00:00').getTime()
    },
    limit: 1000
  }
});
```

#### Get Time Series for Experiment

Retrieves all time series for a specific experiment.

```typescript
// Function: experiments.getExperimentTimeSeries
// Input:
const timeSeries = await client.query("experiments.getExperimentTimeSeries", {
  experimentId: "conv:enhancedExperiments:abc123"
});
```

### Validation Rules

#### Create Validation Rule

Creates a new validation rule for experimental data.

```typescript
// Function: experiments.createValidationRule
// Input:
const ruleId = await client.mutation("experiments.createValidationRule", {
  validationRule: {
    parameterName: "viability_percentage",
    ruleType: "range", // range, pattern, comparison, custom
    parameters: {
      min: 0,
      max: 100,
      inclusiveMin: true,
      inclusiveMax: true
    },
    description: "Viability percentage must be between 0 and 100",
    severity: "error" // error, warning, info
  }
});
```

#### Get Validation Rules

Retrieves validation rules for a specific parameter.

```typescript
// Function: experiments.getValidationRules
// Input:
const rules = await client.query("experiments.getValidationRules", {
  parameterName: "viability_percentage"
});
```

#### Validate Experiment Data

Validates experiment or result data against defined rules.

```typescript
// Function: experiments.validateExperimentData
// Input:
const validationResults = await client.query("experiments.validateExperimentData", {
  experimentId: "conv:enhancedExperiments:abc123"
});
```

## HTTP API

The HTTP API provides RESTful access to the Enhanced Experimental Data System for external applications and systems that don't use Convex directly.

### HTTP Authentication

The HTTP API uses JWT token authentication. Include the token in the `Authorization` header:

```
Authorization: Bearer <your_token>
```

### HTTP Response Format

All API endpoints follow the standardized response format:

```json
{
  "status": "success|error",
  "timestamp": "2025-05-18T14:30:00Z",
  "code": 200,
  "message": "Human-readable message",
  "data": { ... },
  "meta": { ... },
  "pagination": { ... },
  "errors": [ ... ]
}
```

### HTTP Endpoints

The HTTP API mirrors the functionality of the Convex API with RESTful endpoints:

#### Experiments

- `GET /api/v1/experiments` - List/search experiments
- `GET /api/v1/experiments/:id` - Get experiment by ID
- `POST /api/v1/experiments` - Create experiment
- `PUT /api/v1/experiments/:id` - Update experiment
- `DELETE /api/v1/experiments/:id` - Delete experiment
- `PATCH /api/v1/experiments/:id/status` - Update experiment status

#### Experiment Results

- `GET /api/v1/experiments/:experiment_id/results` - Get results for experiment
- `POST /api/v1/experiments/:experiment_id/results` - Add result to experiment
- `PUT /api/v1/results/:id` - Update result
- `DELETE /api/v1/results/:id` - Delete result

#### Protocols

- `GET /api/v1/protocols` - List protocols
- `GET /api/v1/protocols/:id` - Get protocol by ID
- `POST /api/v1/protocols` - Create protocol
- `PUT /api/v1/protocols/:id` - Update protocol
- `POST /api/v1/protocols/:id/versions` - Create new protocol version

#### Tissue Types

- `GET /api/v1/tissue-types` - List tissue types
- `GET /api/v1/tissue-types/:id` - Get tissue type by ID
- `POST /api/v1/tissue-types` - Create tissue type

#### Equipment

- `GET /api/v1/equipment-types` - List equipment types
- `POST /api/v1/equipment-types` - Create equipment type
- `GET /api/v1/equipment` - List equipment
- `POST /api/v1/equipment` - Create equipment
- `POST /api/v1/experiments/:experiment_id/equipment` - Link equipment to experiment
- `GET /api/v1/experiments/:experiment_id/equipment` - Get equipment for experiment

#### Time Series

- `POST /api/v1/experiments/:experiment_id/time-series` - Create time series
- `GET /api/v1/experiments/:experiment_id/time-series` - Get time series for experiment
- `GET /api/v1/time-series/:id` - Get time series with data
- `POST /api/v1/time-series/:id/data-points` - Add data point to time series

## Error Handling

API endpoints return appropriate HTTP status codes and error messages:

```json
{
  "status": "error",
  "timestamp": "2025-05-18T14:30:00Z",
  "code": 404,
  "message": "Resource not found",
  "errors": [
    {
      "type": "NotFoundError",
      "message": "Experiment with ID 'exp-999' not found",
      "context": "GET /experiments/exp-999"
    }
  ]
}
```

Common error codes:

- 400: Bad Request - Invalid input parameters
- 401: Unauthorized - Authentication required
- 403: Forbidden - Insufficient permissions
- 404: Not Found - Resource doesn't exist
- 422: Unprocessable Entity - Validation error
- 500: Internal Server Error - Server-side error

All API operations are validated for:

1. Input validity (format, required fields)
2. Resource existence (referenced entities must exist)
3. Access permissions (public/private access control)
4. Data integrity (foreign key relationships)
5. Business rules (status transitions, etc.)

## Pagination and Filtering

List endpoints support pagination with the following parameters:

- `limit` / `per_page`: Number of items per page (default: 20, max: 100)
- `cursor` / `page`: Pagination cursor or page number
- Multiple filter parameters specific to each endpoint

Pagination information is included in the response:

```json
"pagination": {
  "page": 1,
  "per_page": 20,
  "total_items": 42,
  "total_pages": 3,
  "has_next": true,
  "has_prev": false,
  "cursor": "abc123def456",
  "links": {
    "self": "/api/v1/experiments?page=1&per_page=20",
    "first": "/api/v1/experiments?page=1&per_page=20",
    "last": "/api/v1/experiments?page=3&per_page=20",
    "next": "/api/v1/experiments?page=2&per_page=20"
  }
}
```