# Phase 2.1: API Documentation Implementation for Roo Code

## Objective
Create comprehensive API documentation using OpenAPI/Swagger standards with a focus on micro-tasks for efficient autonomous agent implementation.

## Progress Status
- ✅ Task 1: Install Flask-APISpec and Dependencies - COMPLETE
- ✅ Task 2: Create Core Documentation Infrastructure - COMPLETE
- ✅ Task 5: Implement Interactive Swagger UI - COMPLETE 
- ⏳ Task 3: Document Core API Resources - PARTIALLY COMPLETE
- ❌ Task 4, 6, 7, 8 - NOT STARTED

## Key Files Already Created
- `/api/docs.py` - Core documentation infrastructure
- `/templates/swagger.html` - Swagger UI template
- `/docs/api/openapi.yaml` - Initial OpenAPI specification
- `/api/rdkit_schemas.py` - RDKit-specific schemas

## Current Phase Focus: Resource Schema Documentation

### Task 3A: Molecule Schema Documentation (MICRO-TASK)
```
Agent Task: Create and document the Molecule resource schema.

Inputs:
- Current models in /api/models.py
- OpenAPI specification in /docs/api/openapi.yaml

Steps:
1. Create a MoleculeSchema class in either /api/models.py or a new file /api/schemas.py
2. Define all fields with proper types, validation, and descriptive documentation
3. Add example values for each field
4. Ensure schema is registered with the OpenAPI specification
5. Update spec.path() in /api/docs.py to reference this schema for relevant endpoints

Output:
- Completed Molecule schema with OpenAPI documentation
- Updated path documentation in /api/docs.py for /api/v1/molecules endpoints

Acceptance Criteria:
- Schema includes all Molecule fields (id, name, smiles, inchi, formula, etc.)
- All fields have proper validation, descriptions, and examples
- Schema is used in the OpenAPI path documentation
- Swagger UI shows the schema correctly
```

### Task 3B: Mixture Schema Documentation (MICRO-TASK)
```
Agent Task: Create and document the Mixture resource schema.

Inputs:
- Current models in /api/models.py
- OpenAPI specification in /docs/api/openapi.yaml

Steps:
1. Create a MixtureSchema class in either /api/models.py or a new file /api/schemas.py
2. Define all fields with proper types, validation, and descriptive documentation
3. Create a MixtureComponentSchema for mixture components
4. Add example values for each field
5. Ensure schema is registered with the OpenAPI specification
6. Update spec.path() in /api/docs.py to reference this schema for relevant endpoints

Output:
- Completed Mixture and MixtureComponent schemas with OpenAPI documentation
- Updated path documentation in /api/docs.py for /api/v1/mixtures endpoints

Acceptance Criteria:
- Schema includes all Mixture fields (id, name, description, components, etc.)
- Component relationship is properly modeled
- All fields have proper validation, descriptions, and examples
- Schema is used in the OpenAPI path documentation
- Swagger UI shows the schema correctly
```

### Task 3C: Experiment Schema Documentation (MICRO-TASK)
```
Agent Task: Create and document the Experiment resource schema.

Inputs:
- Current models in /api/models.py
- OpenAPI specification in /docs/api/openapi.yaml

Steps:
1. Create an ExperimentSchema class in either /api/models.py or a new file /api/schemas.py
2. Define all fields with proper types, validation, and descriptive documentation
3. Add example values for each field
4. Ensure schema is registered with the OpenAPI specification
5. Update spec.path() in /api/docs.py to reference this schema for relevant endpoints

Output:
- Completed Experiment schema with OpenAPI documentation
- Updated path documentation in /api/docs.py for /api/v1/experiments endpoints

Acceptance Criteria:
- Schema includes all Experiment fields (id, mixture_id, property_type_id, etc.)
- All fields have proper validation, descriptions, and examples
- Schema is used in the OpenAPI path documentation
- Swagger UI shows the schema correctly
```

### Task 3D: Response Schema Standardization (MICRO-TASK)
```
Agent Task: Create standardized response schemas for API endpoints.

Inputs:
- Current response handling in /api/utils.py and /api/resources.py
- OpenAPI specification in /docs/api/openapi.yaml

Steps:
1. Create or update StandardResponseSchema in /api/models.py or /api/schemas.py
2. Define schemas for success, error, pagination, and metadata responses
3. Add example values for each response type
4. Update spec.path() in /api/docs.py to use these schemas for all endpoints

Output:
- Standardized response schemas with examples
- Updated path documentation using these schemas

Acceptance Criteria:
- Response schemas handle success, error, paginated, and batch responses
- All schemas have proper validation, descriptions, and examples
- Schemas are consistently used across all API paths
- Swagger UI shows the schemas correctly
```

### Task 4A: RDKit Endpoints Documentation (MICRO-TASK)
```
Agent Task: Document RDKit-specific API endpoints using existing schemas.

Inputs:
- Existing RDKit schemas in /api/rdkit_schemas.py
- Current implementation in /api/rdkit_resources.py
- OpenAPI specification in /docs/api/openapi.yaml

Steps:
1. Review and update existing RDKit endpoint documentation in /api/docs.py
2. Ensure all parameters are properly documented
3. Reference the appropriate schemas from rdkit_schemas.py
4. Add detailed descriptions and examples for scientific context

Output:
- Comprehensive documentation for all RDKit endpoints
- Updated path entries in /api/docs.py

Acceptance Criteria:
- All RDKit endpoints are fully documented
- Proper schemas are referenced for request/response bodies
- Scientific context is provided for calculations
- Swagger UI renders the documentation correctly
```

### Task 4B: Predictive Models Documentation (MICRO-TASK)
```
Agent Task: Document predictive model endpoints and schemas.

Inputs:
- Current implementation in /api/predictive_models.py and /api/predictive_models_resources.py
- OpenAPI specification in /docs/api/openapi.yaml

Steps:
1. Create schemas for predictive model requests and responses
2. Document all predictive model endpoints in /api/docs.py
3. Add detailed descriptions of model functionality and parameters
4. Include examples of prediction requests and responses

Output:
- Schemas for predictive model resources
- Comprehensive documentation for prediction endpoints
- Updated path entries in /api/docs.py

Acceptance Criteria:
- All prediction endpoints are fully documented
- Proper schemas for model training and prediction
- Scientific context is provided for models
- Swagger UI renders the documentation correctly
```

## Implementation Instructions for Roo Code PM

1. **Sequential Approach**: Assign these micro-tasks in sequence, starting with Task 3A
2. **Schema First**: Focus on completing all schema tasks before moving to endpoint documentation
3. **Testing**: After each schema is completed, verify it in the Swagger UI before proceeding
4. **Consolidation**: Periodically regenerate the OpenAPI YAML file to ensure consistency
5. **Dependencies**: Note that each task builds on previous ones, especially the schema tasks

## Testing Protocol

After each micro-task:
1. Restart the Flask application 
2. Navigate to the Swagger UI at `/swagger-ui/`
3. Verify that the schema appears correctly
4. Validate that endpoints reference the schema properly
5. Test any example values to ensure they are representative

## Specialized Agent Selection

- Schema Tasks (3A-3D): Assign to Schema Specialist Agent
- RDKit Documentation (4A): Assign to Scientific Documentation Agent
- Predictive Models (4B): Assign to ML Documentation Agent

## Next Phase Planning

After these tasks are complete, use this same micro-task pattern for:
1. Authentication documentation
2. Specific examples in multiple languages
3. Error documentation enhancement