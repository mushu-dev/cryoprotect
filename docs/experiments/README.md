# CryoProtect Experimental Data System

This directory contains documentation and resources for the enhanced experimental data system of CryoProtect.

## Overview

The enhanced experimental data system provides capabilities for:

- Design and execution of cryoprotection experiments
- Protocol generation and versioning
- Result tracking with uncertainty quantification
- Time series data collection and analysis
- Equipment tracking and integration
- Data validation and quality assurance

## System Architecture

The experimental data system is implemented using Convex and follows a modular, layered architecture:

1. **Schema Layer**: Convex schema definitions for experimental data
2. **Type Definitions**: TypeScript interfaces defining the data models
3. **CRUD Operations**: Convex query and mutation functions
4. **Validation Layer**: Input validation and data integrity checks
5. **Helper Functions**: Utilities for data processing and expansion
6. **API Layer**: Accessible APIs for frontend integration

## Components

### Convex Schema

The enhanced schema consists of the following core tables:

- `protocols`: Versioned protocol templates and instances
- `tissueTypes`: Biological samples used in experiments
- `enhancedExperiments`: Main table for comprehensive experiment records
- `enhancedExperimentResults`: Results with uncertainty quantification
- `equipmentTypes`: Categories of laboratory equipment
- `equipment`: Specific equipment instances
- `experimentEquipment`: Links between experiments and equipment
- `timeSeries`: Time-dependent data definitions
- `timeSeriesData`: Individual time-series data points
- `validationRules`: Rules for validating experimental data

### Key Features

#### Protocol Templates and Versioning

The protocol system supports:

- Creating reusable protocol templates
- Detailed step-by-step procedures
- Protocol versioning with parent-child relationships
- Public and private protocols

#### Uncertainty Quantification

Results can include uncertainty information:

- Standard deviation for single-value uncertainty
- Min/max ranges for bounded values
- Confidence intervals with specified confidence levels

#### Time Series Data

Time-dependent data is fully supported:

- Create multiple time series per experiment
- Record individual data points with timestamps
- Include uncertainty in time series measurements
- Metadata for time series context

#### Equipment Tracking

Link experiments to specific equipment:

- Categorize equipment by type
- Track calibration and maintenance information
- Record equipment parameters used in experiments
- Link multiple equipment items to each experiment

### Extended Models

The enhanced system extends or replaces several existing models:

1. `EnhancedExperiment` extends the basic `Experiment` model with:
   - Protocol references
   - Detailed environmental conditions
   - Standardized units
   - Version tracking
   - Provenance information
   - Tagging system

2. `EnhancedExperimentResult` extends the basic `ExperimentResult` model with:
   - Uncertainty quantification
   - Detailed provenance
   - Tissue type references
   - Source documentation

## Usage Examples

### Creating a Protocol Template

```typescript
// Create a new protocol template
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
        temperatureUnit: "C",
      }
    ],
    isTemplate: true,
    category: "vitrification",
    public: true
  }
});
```

### Creating an Experiment with Protocol

```typescript
// Create a new experiment using a protocol template
const experimentId = await client.mutation("experiments.createEnhancedExperiment", {
  experiment: {
    name: "DMSO Vitrification Test #5",
    description: "Testing vitrification with DMSO concentration series",
    protocolId,  // Reference to the protocol template
    mixtureId,   // Reference to the DMSO mixture
    temperature: -196,
    temperatureUnit: "C",
    coolingRate: -200,
    coolingRateUnit: "C/min",
    thawingRate: 10,
    thawingRateUnit: "C/min",
    parameters: {
      sampleCount: 10,
      controlGroup: true
    },
    status: "in-progress",
    tags: ["vitrification", "dmso", "mouse-embryos"],
    date: Date.now()
  }
});
```

### Recording Results with Uncertainty

```typescript
// Add a result with uncertainty quantification
const resultId = await client.mutation("experiments.createEnhancedExperimentResult", {
  result: {
    experimentId,
    tissueTypeId,
    parameterName: "viability_percentage",
    value: 87.5,
    numericValue: 87.5,
    units: "%",
    uncertainty: {
      type: "standard_deviation",
      value: 4.2
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

### Working with Time Series Data

```typescript
// Create a time series for temperature measurements
const timeSeriesId = await client.mutation("experiments.createTimeSeries", {
  experimentId,
  name: "Cooling Curve",
  parameterName: "temperature",
  units: "C"
});

// Add time series data points
await client.mutation("experiments.addTimeSeriesDataPoint", {
  timeSeriesId,
  timestamp: Date.now(),
  value: -5.2,
  uncertainty: 0.1
});
```

### Querying Complete Experiment Data

```typescript
// Get a complete experiment with all related data
const experiment = await client.query("experiments.getEnhancedExperiment", {
  experimentId,
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

## Development Guide

### Implementation Structure

The enhanced experimental system is implemented in the following files:

```
convex/
├── schema/
│   ├── convex_schema.ts         # Main schema
│   └── enhanced_experiment_schema.ts # Enhanced schema
└── experiments/
    ├── types.ts                 # Original experiment types
    ├── enhanced_types.ts        # Enhanced experiment types
    ├── enhanced_validation.ts   # Input validation
    ├── enhanced_helpers.ts      # Helper functions
    ├── enhanced_experiments.ts  # Core experiment operations
    ├── enhanced_experiment_results.ts # Result operations
    ├── protocols.ts             # Protocol operations
    └── enhanced_index.ts        # Unified exports
```

### Adding New Functionality

To extend the experimental data system:

1. Define new types in `enhanced_types.ts`
2. Add validation functions in `enhanced_validation.ts`
3. Implement helper functions in `enhanced_helpers.ts`
4. Create CRUD operations in a new or existing module
5. Export from `enhanced_index.ts`
6. Update the schema if necessary

### Working with Time Series Data

The time series system provides flexible support for temporal data:

1. Create a time series definition with `createTimeSeries`
2. Add data points with `addTimeSeriesDataPoint`
3. Query time series data with `getTimeSeriesWithData`
4. Analyze time series data with specialized methods

### Equipment Integration

To integrate laboratory equipment:

1. Define equipment types with `createEquipmentType`
2. Register specific equipment with `createEquipment`
3. Link equipment to experiments with `createExperimentEquipment`
4. Query equipment used in experiments

## Integration with Frontend

To integrate with the frontend:

1. Import experimental data types:
   ```typescript
   import { 
     EnhancedExperiment,
     EnhancedExperimentResult,
     Protocol,
     TissueType
   } from '../convex/experiments/enhanced_index';
   ```

2. Use the Convex hooks for data access:
   ```typescript
   import { useQuery, useMutation } from "convex/react";
   import { api } from "../convex/_generated/api";
   
   // Query experiments
   const experiments = useQuery(api.experiments.listEnhancedExperiments);
   
   // Create experiment
   const createExperiment = useMutation(api.experiments.createEnhancedExperiment);
   ```

## Resources

- [Experimental Data Enhancement System](../../EXPERIMENTAL_DATA_ENHANCEMENT_SYSTEM.md)
- [Schema Design](./schema_design.md)
- [API Reference](./api.md)
- [Scientific Methods](./scientific_methods.md)