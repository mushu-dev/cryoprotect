"""
CryoProtect Analyzer API - Schema Definitions

This module contains Marshmallow schemas for API validation and documentation.
These schemas define the structure, validation rules, and documentation for API
request and response objects.
"""

from marshmallow import Schema, fields, validate, validates, ValidationError
import re
from typing import List, Dict, Any, Optional

from api.models import FlexibleDateTime

# Lab verification fields
lab_verification_fields = {
    'id': fields.String,
    'experiment_id': fields.String,
    'verification_status': fields.String,
    'verifier': fields.String,
    'equipment_used': fields.String,
    'comments': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601'),
    'updated_at': FlexibleDateTime(dt_format='iso8601')
}

verification_stats_fields = {
    'total_count': fields.Integer,
    'verified_count': fields.Integer,
    'pending_count': fields.Integer,
    'rejected_count': fields.Integer,
    'verification_rate': fields.Float,
    'by_equipment': fields.Raw,
    'by_verifier': fields.Raw
}

# Lab verification schemas
class LabVerificationSchema(Schema):
    """
    Schema for lab verification data.
    
    This schema defines the structure for validation data associated with
    experiments, tracking verification status and metadata.
    """
    id = fields.UUID(
        description="Unique identifier for the verification record",
        example="623e4567-e89b-12d3-a456-426614174005"
    )
    experiment_id = fields.UUID(
        required=True,
        description="ID of the experiment being verified",
        example="523e4567-e89b-12d3-a456-426614174004"
    )
    verification_status = fields.String(
        required=True,
        description="Status of verification",
        example="verified",
        validate=validate.OneOf(["pending", "verified", "rejected"])
    )
    verifier = fields.String(
        required=True,
        description="User ID or name of the person performing verification",
        example="a23e4567-e89b-12d3-a456-426614174000"
    )

# User profile schema
class UserProfileSchema(Schema):
    """
    Schema for user profile creation and update.
    """
    email = fields.Email(
        required=True,
        description="User email address",
        example="user@example.com"
    )
    name = fields.String(
        required=False,
        description="User's full name",
        example="Jane Doe"
    )

# User profile schema
class UserProfileSchema(Schema):
    """
    Schema for user profile creation and update.
    """
    email = fields.Email(
        required=True,
        description="User email address",
        example="user@example.com"
    )
    name = fields.String(
        required=False,
        description="User's full name",
        example="Jane Doe"
    )
# Property comparison request schema
class PropertyComparisonRequestSchema(Schema):
    """
    Schema for property comparison request.
    """
    ids = fields.List(
        fields.String(),
        required=True,
        description="List of entity IDs (molecules or mixtures) to compare"
    )
    comments = fields.String(
        description="Additional notes or comments about verification",
        example="Results match prediction within 2% margin of error"
    )
    created_at = fields.DateTime(
        description="Timestamp of when the verification was recorded",
        example="2023-04-05T14:20:00Z"
    )
    updated_at = fields.DateTime(
        description="Timestamp of when the verification was last updated",
        example="2023-04-06T09:15:00Z"
    )

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