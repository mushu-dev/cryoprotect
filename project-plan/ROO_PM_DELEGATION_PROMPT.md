# Roo PM Delegation Prompt

```
Agent = ProjectManager

I need your help implementing the CryoProtect v2 project according to our Direct Execution Guide. This document contains precise implementation instructions that remove all guesswork from the process.

## Instructions

1. Review the `/project-plan/DIRECT_EXECUTION_GUIDE.md` document thoroughly to understand the current API documentation implementation tasks
2. For each task in the guide:
   - Assign it to the appropriate specialized agent (Schema Specialist, Documentation Engineer, etc.)
   - Provide the exact file paths and code from the guide in your delegation
   - Specify the success criteria exactly as defined in the Expected Outcome section
   - Request verification that matches the verification steps in the guide
3. Monitor the implementation for precise compliance with the specification
4. Verify completed work against the Expected Outcome criteria before proceeding to the next task
5. Maintain a running status report of completed and pending tasks

## Task Delegation Format

When delegating a task, use this format:

```
TASK ASSIGNMENT: [Task Name]

IMPLEMENTING AGENT: [Agent Type]

REFERENCE: Direct Execution Guide - Task [Number]

FILES TO MODIFY:
- [Exact file path]

IMPLEMENTATION INSTRUCTIONS:
[Paste the exact code and instructions from the guide]

VERIFICATION STEPS:
1. [Steps from the guide's verification section]
2. ...

EXPECTED OUTCOME:
[Copy from the guide's Expected Outcome section]
```

## Current Implementation Focus

The current focus is on API documentation tasks in Phase 2.1:
- TASK_1: Create Core Resource Schemas (MoleculeSchema)
- TASK_2: Create Mixture Schema (MixtureSchema, MixtureComponentSchema)
- TASK_3: Create Experiment Schema (ExperimentSchema, PropertyTypeSchema)
- TASK_4: Create Standardized Response Schemas (ErrorResponseSchema, SuccessResponseSchema)

Once these schemas are implemented, we'll proceed to the documentation integration tasks.

Please begin by reviewing the guide and delegating TASK_1 to a Schema Specialist Agent using the exact implementation instructions.
```
TASK ASSIGNMENT: Create Experiment Schema

IMPLEMENTING AGENT: Schema Specialist

REFERENCE: Direct Execution Guide - Task 3

FILES TO MODIFY:
- /api/schemas.py
- /api/docs.py

IMPLEMENTATION INSTRUCTIONS:
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

VERIFICATION STEPS:
1. Confirm that `/api/schemas.py` contains the new PropertyTypeSchema, ExperimentSchema, ExperimentListSchema, and ExperimentCreateSchema classes as specified.
2. Confirm that `/api/docs.py` imports the new schemas as shown.
3. Confirm that `/api/docs.py` includes the experiment endpoint documentation block exactly as specified.

EXPECTED OUTCOME:
- Updated `/api/schemas.py` with Experiment schemas
- Updated imports in `/api/docs.py`
- Added experiment endpoints documentation to `/api/docs.py`
## Status Report Update (2025-04-21 17:29 MDT)
- TASK_1: Complete
- TASK_2: Complete
- TASK_3: Complete (Experiment and PropertyType schemas implemented, docs.py updated per Direct Execution Guide)
- TASK_4: Pending