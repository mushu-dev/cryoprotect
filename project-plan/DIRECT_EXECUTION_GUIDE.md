# CryoProtect v2 Direct Execution Guide

## Purpose
This document provides explicit, step-by-step instructions for Roo agents to implement specific components of the CryoProtect v2 project without any guesswork. Each implementation task includes exact file locations, code patterns, and expected outcomes.

## How to Use This Guide

As the Project Manager, use the following prompt pattern with Roo:

```
I want you to implement TASK_NUMBER from the Direct Execution Guide. Please follow the exact implementation instructions without deviation. Focus only on writing the specified code in the specified files.
```

## Current Implementation Focus: API Documentation Phase

### TASK_1: Create Core Resource Schemas

**Objective**: Create Marshmallow schemas for the Molecule resource in a dedicated schema file.

**Files to Create/Modify**:
1. Create file: `/api/schemas.py`

**Exact Implementation**:

1. Create `/api/schemas.py` with the following content:

```python
"""
CryoProtect Analyzer API - Schema Definitions

This module contains Marshmallow schemas for API validation and documentation.
These schemas define the structure, validation rules, and documentation for API
request and response objects.
"""

from marshmallow import Schema, fields, validate, validates, ValidationError
import re
from typing import List, Dict, Any, Optional

class MoleculeSchema(Schema):
    """
    Schema for Molecule resource.
    
    This schema defines the structure and validation for molecule data used in
    API requests and responses. It includes fields for molecular identifiers,
    structure data, and metadata.
    """
    id = fields.UUID(
        description="Unique identifier for the molecule",
        example="123e4567-e89b-12d3-a456-426614174000"
    )
    name = fields.String(
        required=True,
        description="Name of the molecule",
        example="Glycerol",
        validate=validate.Length(min=1, max=255)
    )
    smiles = fields.String(
        required=True,
        description="SMILES representation of the molecule",
        example="C(C(CO)O)O"
    )
    inchi = fields.String(
        description="InChI identifier for the molecule",
        example="InChI=1S/C3H8O3/c4-1-3(6)2-5/h3,5-6H,1-2H2,4H"
    )
    inchi_key = fields.String(
        description="InChI key for the molecule",
        example="PEDCQBHIVMGVHV-UHFFFAOYSA-N"
    )
    molecular_formula = fields.String(
        description="Molecular formula",
        example="C3H8O3"
    )
    description = fields.String(
        description="Additional information about the molecule",
        example="Common cryoprotectant used in various freezing protocols"
    )
    tags = fields.List(
        fields.String(),
        description="List of tags or categories for the molecule",
        example=["cryoprotectant", "polyol", "FDA-approved"]
    )
    created_at = fields.DateTime(
        description="Timestamp of when the molecule was created",
        example="2023-04-01T12:00:00Z"
    )
    updated_at = fields.DateTime(
        description="Timestamp of when the molecule was last updated",
        example="2023-04-15T14:30:00Z"
    )
    user_id = fields.UUID(
        description="ID of the user who created the molecule",
        example="a23e4567-e89b-12d3-a456-426614174000"
    )
    
    @validates("smiles")
    def validate_smiles(self, value):
        """Validate SMILES string."""
        if not value or not isinstance(value, str):
            raise ValidationError("SMILES must be a non-empty string")
        
        # Basic validation - could be enhanced with RDKit validation
        if not re.match(r'^[A-Za-z0-9@+\-\[\]\(\)\\\/%\.=#$:,~]+$', value):
            raise ValidationError("SMILES contains invalid characters")
        
        return True

class MoleculeListSchema(Schema):
    """Schema for a list of molecules with pagination."""
    molecules = fields.List(
        fields.Nested(MoleculeSchema),
        required=True,
        description="List of molecule objects"
    )
    total = fields.Integer(
        required=True,
        description="Total number of molecules matching the query",
        example=42
    )
    page = fields.Integer(
        required=True,
        description="Current page number",
        example=1
    )
    per_page = fields.Integer(
        required=True,
        description="Number of items per page",
        example=20
    )
    
class MoleculeCreateSchema(Schema):
    """Schema for creating a new molecule."""
    name = fields.String(
        required=True,
        description="Name of the molecule",
        example="Glycerol",
        validate=validate.Length(min=1, max=255)
    )
    smiles = fields.String(
        required=True,
        description="SMILES representation of the molecule",
        example="C(C(CO)O)O"
    )
    inchi = fields.String(
        description="InChI identifier for the molecule",
        example="InChI=1S/C3H8O3/c4-1-3(6)2-5/h3,5-6H,1-2H2,4H"
    )
    molecular_formula = fields.String(
        description="Molecular formula",
        example="C3H8O3"
    )
    description = fields.String(
        description="Additional information about the molecule",
        example="Common cryoprotectant used in various freezing protocols"
    )
    tags = fields.List(
        fields.String(),
        description="List of tags or categories for the molecule",
        example=["cryoprotectant", "polyol", "FDA-approved"]
    )
```

2. Modify `/api/docs.py` to reference the new schemas by adding these imports at the top of the file:

```python
from api.schemas import MoleculeSchema, MoleculeListSchema, MoleculeCreateSchema
```

3. Update the molecule endpoint documentation in `/api/docs.py` to reference the schemas:

```python
# Molecule endpoints
spec.path(
    path="/api/v1/molecules",
    operations={
        'get': {
            'summary': 'Get a list of molecules',
            'description': 'Returns a paginated list of molecules',
            'parameters': [
                {'name': 'page', 'in': 'query', 'schema': {'type': 'integer'}, 'description': 'Page number'},
                {'name': 'per_page', 'in': 'query', 'schema': {'type': 'integer'}, 'description': 'Items per page'},
                {'name': 'search', 'in': 'query', 'schema': {'type': 'string'}, 'description': 'Search term'}
            ],
            'responses': {
                '200': {
                    'description': 'List of molecules',
                    'content': {
                        'application/json': {
                            'schema': MoleculeListSchema
                        }
                    }
                }
            }
        },
        'post': {
            'summary': 'Create a new molecule',
            'description': 'Creates a new molecule from the provided data',
            'requestBody': {
                'required': True,
                'content': {
                    'application/json': {
                        'schema': MoleculeCreateSchema
                    }
                }
            },
            'responses': {
                '201': {
                    'description': 'Molecule created successfully',
                    'content': {
                        'application/json': {
                            'schema': MoleculeSchema
                        }
                    }
                },
                '400': {
                    'description': 'Invalid request data'
                }
            }
        }
    }
)
```

**Expected Outcome**:
- New file `/api/schemas.py` with Molecule schemas
- Updated imports in `/api/docs.py`
- Updated molecule endpoints documentation in `/api/docs.py`

### TASK_2: Create Mixture Schema

**Objective**: Add Mixture schemas to the schema definitions.

**Files to Modify**:
1. Update file: `/api/schemas.py`

**Exact Implementation**:

Add the following code to the end of `/api/schemas.py`:

```python
class MixtureComponentSchema(Schema):
    """
    Schema for a mixture component.
    
    This schema defines the structure and validation for mixture components,
    representing the molecules and their concentrations within a mixture.
    """
    id = fields.UUID(
        description="Unique identifier for the mixture component",
        example="223e4567-e89b-12d3-a456-426614174001"
    )
    mixture_id = fields.UUID(
        required=True,
        description="ID of the parent mixture",
        example="323e4567-e89b-12d3-a456-426614174002"
    )
    molecule_id = fields.UUID(
        required=True,
        description="ID of the molecule in this component",
        example="123e4567-e89b-12d3-a456-426614174000"
    )
    concentration = fields.Float(
        required=True,
        description="Concentration of the molecule in percentage or molarity",
        example=10.5,
        validate=validate.Range(min=0)
    )
    concentration_unit = fields.String(
        required=True,
        description="Unit of concentration (e.g., percentage, molarity)",
        example="percentage",
        validate=validate.OneOf(["percentage", "molarity", "mg/ml", "mM", "µM"])
    )
    molecule = fields.Nested(
        MoleculeSchema,
        description="The molecule data (included in GET responses)",
        dump_only=True
    )

class MixtureSchema(Schema):
    """
    Schema for Mixture resource.
    
    This schema defines the structure and validation for mixture data used in
    API requests and responses. It includes fields for mixture metadata and
    components.
    """
    id = fields.UUID(
        description="Unique identifier for the mixture",
        example="323e4567-e89b-12d3-a456-426614174002"
    )
    name = fields.String(
        required=True,
        description="Name of the mixture",
        example="Glycerol-DMSO Mixture",
        validate=validate.Length(min=1, max=255)
    )
    description = fields.String(
        description="Description of the mixture and its uses",
        example="Standard cryopreservation mixture for cell freezing"
    )
    components = fields.List(
        fields.Nested(MixtureComponentSchema),
        description="List of components in the mixture",
        dump_only=True
    )
    tags = fields.List(
        fields.String(),
        description="List of tags or categories for the mixture",
        example=["cell freezing", "10% glycerol", "rapid cooling"]
    )
    created_at = fields.DateTime(
        description="Timestamp of when the mixture was created",
        example="2023-04-02T14:30:00Z"
    )
    updated_at = fields.DateTime(
        description="Timestamp of when the mixture was last updated",
        example="2023-04-16T09:45:00Z"
    )
    user_id = fields.UUID(
        description="ID of the user who created the mixture",
        example="a23e4567-e89b-12d3-a456-426614174000"
    )

class MixtureListSchema(Schema):
    """Schema for a list of mixtures with pagination."""
    mixtures = fields.List(
        fields.Nested(MixtureSchema),
        required=True,
        description="List of mixture objects"
    )
    total = fields.Integer(
        required=True,
        description="Total number of mixtures matching the query",
        example=15
    )
    page = fields.Integer(
        required=True,
        description="Current page number",
        example=1
    )
    per_page = fields.Integer(
        required=True,
        description="Number of items per page",
        example=20
    )

class MixtureCreateSchema(Schema):
    """Schema for creating a new mixture."""
    name = fields.String(
        required=True,
        description="Name of the mixture",
        example="Glycerol-DMSO Mixture",
        validate=validate.Length(min=1, max=255)
    )
    description = fields.String(
        description="Description of the mixture and its uses",
        example="Standard cryopreservation mixture for cell freezing"
    )
    components = fields.List(
        fields.Dict(keys=fields.String(), values=fields.Raw()),
        required=True,
        description="List of component objects with molecule_id, concentration, and concentration_unit",
        example=[
            {
                "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
                "concentration": 10.0,
                "concentration_unit": "percentage"
            },
            {
                "molecule_id": "223e4567-e89b-12d3-a456-426614174001",
                "concentration": 5.0,
                "concentration_unit": "percentage"
            }
        ]
    )
    tags = fields.List(
        fields.String(),
        description="List of tags or categories for the mixture",
        example=["cell freezing", "10% glycerol", "rapid cooling"]
    )
```

Update the imports in `/api/docs.py` to include the new schemas:

```python
from api.schemas import (
    MoleculeSchema, MoleculeListSchema, MoleculeCreateSchema,
    MixtureSchema, MixtureListSchema, MixtureCreateSchema
)
```

Add mixture endpoint documentation to `/api/docs.py`:

```python
# Mixture endpoints
spec.path(
    path="/api/v1/mixtures",
    operations={
        'get': {
            'summary': 'Get a list of mixtures',
            'description': 'Returns a paginated list of mixtures',
            'parameters': [
                {'name': 'page', 'in': 'query', 'schema': {'type': 'integer'}, 'description': 'Page number'},
                {'name': 'per_page', 'in': 'query', 'schema': {'type': 'integer'}, 'description': 'Items per page'},
                {'name': 'search', 'in': 'query', 'schema': {'type': 'string'}, 'description': 'Search term'}
            ],
            'responses': {
                '200': {
                    'description': 'List of mixtures',
                    'content': {
                        'application/json': {
                            'schema': MixtureListSchema
                        }
                    }
                }
            }
        },
        'post': {
            'summary': 'Create a new mixture',
            'description': 'Creates a new mixture from the provided data',
            'requestBody': {
                'required': True,
                'content': {
                    'application/json': {
                        'schema': MixtureCreateSchema
                    }
                }
            },
            'responses': {
                '201': {
                    'description': 'Mixture created successfully',
                    'content': {
                        'application/json': {
                            'schema': MixtureSchema
                        }
                    }
                },
                '400': {
                    'description': 'Invalid request data'
                }
            }
        }
    }
)
```

**Expected Outcome**:
- Updated `/api/schemas.py` with Mixture schemas
- Updated imports in `/api/docs.py`
- Added mixture endpoints documentation to `/api/docs.py`

### TASK_3: Create Experiment Schema

**Objective**: Add Experiment schemas to the schema definitions.

**Files to Modify**:
1. Update file: `/api/schemas.py`

**Exact Implementation**:

Add the following code to the end of `/api/schemas.py`:

```python
class PropertyTypeSchema(Schema):
    """
    Schema for property types.
    
    This schema defines the structure for property types that can be measured
    in experiments, such as freezing point, viability, or viscosity.
    """
    id = fields.UUID(
        description="Unique identifier for the property type",
        example="423e4567-e89b-12d3-a456-426614174003"
    )
    name = fields.String(
        required=True,
        description="Name of the property",
        example="Cell Viability",
        validate=validate.Length(min=1, max=255)
    )
    description = fields.String(
        description="Description of the property",
        example="Percentage of cells surviving after freeze/thaw cycle"
    )
    data_type = fields.String(
        required=True,
        description="Data type for this property",
        example="numeric",
        validate=validate.OneOf(["numeric", "text", "boolean"])
    )
    unit = fields.String(
        description="Unit of measurement for numeric properties",
        example="percentage"
    )

class ExperimentSchema(Schema):
    """
    Schema for Experiment resource.
    
    This schema defines the structure and validation for experiment data,
    representing measurements of mixture properties in controlled conditions.
    """
    id = fields.UUID(
        description="Unique identifier for the experiment",
        example="523e4567-e89b-12d3-a456-426614174004"
    )
    mixture_id = fields.UUID(
        required=True,
        description="ID of the mixture used in the experiment",
        example="323e4567-e89b-12d3-a456-426614174002"
    )
    property_type_id = fields.UUID(
        required=True,
        description="ID of the property being measured",
        example="423e4567-e89b-12d3-a456-426614174003"
    )
    numeric_value = fields.Float(
        description="Numeric result of the experiment",
        example=85.2
    )
    text_value = fields.String(
        description="Text result of the experiment",
        example="High viability observed"
    )
    boolean_value = fields.Boolean(
        description="Boolean result of the experiment",
        example=True
    )
    temperature = fields.Float(
        description="Temperature at which the experiment was conducted (°C)",
        example=-80.0
    )
    conditions = fields.Dict(
        keys=fields.String(),
        values=fields.Raw(),
        description="Additional experimental conditions as key-value pairs",
        example={
            "cooling_rate": "1°C/min",
            "storage_time": "48 hours",
            "thawing_method": "rapid water bath"
        }
    )
    notes = fields.String(
        description="Additional notes about the experiment",
        example="Cells showed excellent recovery after 48-hour storage"
    )
    conducted_at = fields.DateTime(
        description="Timestamp of when the experiment was conducted",
        example="2023-04-03T10:15:00Z"
    )
    created_at = fields.DateTime(
        description="Timestamp of when the experiment record was created",
        example="2023-04-03T11:30:00Z"
    )
    user_id = fields.UUID(
        description="ID of the user who conducted the experiment",
        example="a23e4567-e89b-12d3-a456-426614174000"
    )
    property_type = fields.Nested(
        PropertyTypeSchema,
        description="The property type being measured",
        dump_only=True
    )
    mixture = fields.Nested(
        MixtureSchema,
        description="The mixture used in the experiment",
        dump_only=True
    )

class ExperimentListSchema(Schema):
    """Schema for a list of experiments with pagination."""
    experiments = fields.List(
        fields.Nested(ExperimentSchema),
        required=True,
        description="List of experiment objects"
    )
    total = fields.Integer(
        required=True,
        description="Total number of experiments matching the query",
        example=28
    )
    page = fields.Integer(
        required=True,
        description="Current page number",
        example=1
    )
    per_page = fields.Integer(
        required=True,
        description="Number of items per page",
        example=20
    )

class ExperimentCreateSchema(Schema):
    """Schema for creating a new experiment."""
    mixture_id = fields.UUID(
        required=True,
        description="ID of the mixture used in the experiment",
        example="323e4567-e89b-12d3-a456-426614174002"
    )
    property_type_id = fields.UUID(
        required=True,
        description="ID of the property being measured",
        example="423e4567-e89b-12d3-a456-426614174003"
    )
    numeric_value = fields.Float(
        description="Numeric result of the experiment",
        example=85.2
    )
    text_value = fields.String(
        description="Text result of the experiment",
        example="High viability observed"
    )
    boolean_value = fields.Boolean(
        description="Boolean result of the experiment",
        example=True
    )
    temperature = fields.Float(
        description="Temperature at which the experiment was conducted (°C)",
        example=-80.0
    )
    conditions = fields.Dict(
        keys=fields.String(),
        values=fields.Raw(),
        description="Additional experimental conditions as key-value pairs"
    )
    notes = fields.String(
        description="Additional notes about the experiment"
    )
    conducted_at = fields.DateTime(
        description="Timestamp of when the experiment was conducted"
    )
```

Update the imports in `/api/docs.py` again:

```python
from api.schemas import (
    MoleculeSchema, MoleculeListSchema, MoleculeCreateSchema,
    MixtureSchema, MixtureListSchema, MixtureCreateSchema,
    ExperimentSchema, ExperimentListSchema, ExperimentCreateSchema
)
```

Add experiment endpoint documentation to `/api/docs.py`:

```python
# Experiment endpoints
spec.path(
    path="/api/v1/mixtures/{mixture_id}/experiments",
    operations={
        'get': {
            'summary': 'Get experiments for a mixture',
            'description': 'Returns a list of experiments for a mixture',
            'parameters': [
                {'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string', 'format': 'uuid'}}
            ],
            'responses': {
                '200': {
                    'description': 'List of experiments',
                    'content': {
                        'application/json': {
                            'schema': ExperimentListSchema
                        }
                    }
                }
            }
        },
        'post': {
            'summary': 'Record an experiment for a mixture',
            'description': 'Record a new experiment for a mixture',
            'parameters': [
                {'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string', 'format': 'uuid'}}
            ],
            'requestBody': {
                'required': True,
                'content': {
                    'application/json': {
                        'schema': ExperimentCreateSchema
                    }
                }
            },
            'responses': {
                '201': {
                    'description': 'Experiment recorded successfully',
                    'content': {
                        'application/json': {
                            'schema': ExperimentSchema
                        }
                    }
                },
                '400': {
                    'description': 'Invalid request data'
                }
            }
        }
    }
)
```

**Expected Outcome**:
- Updated `/api/schemas.py` with Experiment schemas
- Updated imports in `/api/docs.py`
- Added experiment endpoints documentation to `/api/docs.py`

### TASK_4: Create Standardized Response Schemas

**Objective**: Create standardized response schemas for consistent API responses.

**Files to Modify**:
1. Update file: `/api/schemas.py`
2. Update file: `/api/docs.py`

**Exact Implementation**:

Add the following code at the beginning of `/api/schemas.py` (after the imports):

```python
class ErrorResponseSchema(Schema):
    """
    Schema for standardized error responses.
    
    This schema defines the structure of error responses returned by the API,
    ensuring consistent error reporting across all endpoints.
    """
    status = fields.String(
        required=True,
        description="Status indicator, always 'error' for error responses",
        example="error"
    )
    message = fields.String(
        required=True,
        description="Human-readable error message",
        example="The requested resource was not found"
    )
    code = fields.String(
        description="Error code for programmatic identification",
        example="RESOURCE_NOT_FOUND"
    )
    details = fields.Raw(
        description="Additional error details, varies by error type",
        example={
            "resource_type": "molecule",
            "resource_id": "123e4567-e89b-12d3-a456-426614174000"
        }
    )
    timestamp = fields.DateTime(
        required=True,
        description="Timestamp of when the error occurred",
        example="2023-04-20T14:32:15Z"
    )

class MetadataSchema(Schema):
    """
    Schema for response metadata.
    
    This schema defines standard metadata included with API responses,
    such as pagination information, request processing time, and API version.
    """
    pagination = fields.Dict(
        keys=fields.String(),
        values=fields.Raw(),
        description="Pagination information for list responses",
        example={
            "page": 1,
            "per_page": 20,
            "total": 42,
            "pages": 3
        }
    )
    processing_time = fields.Float(
        description="Request processing time in milliseconds",
        example=45.2
    )
    api_version = fields.String(
        description="API version that processed the request",
        example="1.0.0"
    )
    request_id = fields.String(
        description="Unique identifier for the request (for troubleshooting)",
        example="req_7f8d9e6b-c987-4321-ba98-76c5d4e3f2b1"
    )

class SuccessResponseSchema(Schema):
    """
    Schema for standardized success responses.
    
    This schema defines the structure of success responses returned by the API,
    ensuring consistent response formatting across all endpoints.
    """
    status = fields.String(
        required=True,
        description="Status indicator, always 'success' for successful responses",
        example="success"
    )
    data = fields.Raw(
        description="Response data, varies by endpoint",
        example={"id": "123e4567-e89b-12d3-a456-426614174000", "name": "Glycerol"}
    )
    message = fields.String(
        description="Optional success message",
        example="Resource created successfully"
    )
    metadata = fields.Nested(
        MetadataSchema,
        description="Additional metadata about the response"
    )
    timestamp = fields.DateTime(
        required=True,
        description="Timestamp of the response",
        example="2023-04-20T14:33:22Z"
    )
```

Add this function to the end of `/api/docs.py`:

```python
def add_standard_responses(operations):
    """
    Add standard responses to API operations.
    
    Args:
        operations: Dictionary of operations to update
    
    Returns:
        Updated operations dictionary with standard responses
    """
    for operation in operations.values():
        if 'responses' not in operation:
            operation['responses'] = {}
            
        # Add 400 Bad Request if not present
        if '400' not in operation['responses']:
            operation['responses']['400'] = {
                'description': 'Bad Request - Invalid input data',
                'content': {
                    'application/json': {
                        'schema': ErrorResponseSchema
                    }
                }
            }
            
        # Add 401 Unauthorized if not present
        if '401' not in operation['responses']:
            operation['responses']['401'] = {
                'description': 'Unauthorized - Authentication required',
                'content': {
                    'application/json': {
                        'schema': ErrorResponseSchema
                    }
                }
            }
            
        # Add 404 Not Found for single resource operations
        if '{' in operation.get('path', ''):
            if '404' not in operation['responses']:
                operation['responses']['404'] = {
                    'description': 'Not Found - Resource does not exist',
                    'content': {
                        'application/json': {
                            'schema': ErrorResponseSchema
                        }
                    }
                }
                
        # Add 500 Server Error if not present
        if '500' not in operation['responses']:
            operation['responses']['500'] = {
                'description': 'Internal Server Error',
                'content': {
                    'application/json': {
                        'schema': ErrorResponseSchema
                    }
                }
            }
    
    return operations
```

Update the imports in `/api/docs.py` to include the standard response schemas:

```python
from api.schemas import (
    ErrorResponseSchema, SuccessResponseSchema, MetadataSchema,
    MoleculeSchema, MoleculeListSchema, MoleculeCreateSchema,
    MixtureSchema, MixtureListSchema, MixtureCreateSchema,
    ExperimentSchema, ExperimentListSchema, ExperimentCreateSchema
)
```

Modify all `spec.path()` calls in `/api/docs.py` to use the `add_standard_responses` function. For example:

```python
# Molecule endpoints
spec.path(
    path="/api/v1/molecules",
    operations=add_standard_responses({
        'get': {
            'summary': 'Get a list of molecules',
            # rest of the code remains the same
        },
        'post': {
            'summary': 'Create a new molecule',
            # rest of the code remains the same
        }
    })
)
```

Repeat this pattern for all other endpoint definitions.

**Expected Outcome**:
- Updated `/api/schemas.py` with standardized response schemas
- Updated imports in `/api/docs.py`
- Added utility function for standardized responses
- All API endpoint documentation using standardized error responses

## Next Implementation Phase: Protocol Designer

After completing the API Documentation tasks, the next implementation phase will focus on the Protocol Designer functionality. Specific tasks for that phase will be added here once the API Documentation phase is complete.