"""
CryoProtect Analyzer API - Data Models

This module contains data models for the API, including database models for interacting with Supabase
and request/response schemas for API validation and serialization.

The module is organized into several sections:
- Response fields: Flask-RESTful field definitions for API responses
- Validation schemas: Marshmallow schemas for request validation
- Base models: Generic database operation classes
- Entity models: Specific models for each entity type (molecules, mixtures, etc.)

All schemas include proper validation rules and documentation for OpenAPI generation.
"""

import uuid
from datetime import datetime, date
from typing import Dict, List, Optional, Union, Any, TypeVar, Generic, Type
from flask_restful import fields
from marshmallow import Schema, fields as ma_fields, validate, ValidationError, post_load
from postgrest.exceptions import APIError

# Try to import apispec for schema documentation
try:
    from apispec import APISpec
    from apispec.ext.marshmallow import MarshmallowPlugin
    APISPEC_AVAILABLE = True
except ImportError:
    APISPEC_AVAILABLE = False

from api.utils import get_supabase_client, get_user_id, handle_supabase_error, handle_error

# Type variable for generic model methods
T = TypeVar('T')


class FlexibleDateTime(fields.Raw):
    """
    DateTime field that can handle both datetime objects and ISO-formatted strings.
    
    This field formatter provides flexible handling of datetime values, accepting
    both datetime objects and ISO-formatted strings. It ensures consistent
    datetime formatting in API responses.
    
    Attributes:
        dt_format (str): Format to use for datetime conversion (default: 'iso8601')
    """
    
    def __init__(self, dt_format='iso8601', **kwargs):
        """
        Initialize the FlexibleDateTime field.
        
        Args:
            dt_format (str): Format to use for datetime conversion
            **kwargs: Additional arguments passed to the parent class
        """
        self.dt_format = dt_format
        super(FlexibleDateTime, self).__init__(**kwargs)
    
    def format(self, value):
        """
        Format the datetime value.
        
        Args:
            value: The datetime value to format (datetime object or string)
            
        Returns:
            str: Formatted datetime string or None if input is None
        """
        if value is None:
            return None
        
        # If it's already a string, check if it's ISO format and return as is
        if isinstance(value, str):
            try:
                # Validate it's a proper ISO format by parsing it
                # Replace 'Z' with '+00:00' for compatibility with fromisoformat
                datetime.fromisoformat(value.replace('Z', '+00:00'))
                return value
            except ValueError:
                # If not a valid ISO format, return as is
                return value
        
        # If it's a datetime, format it
        if isinstance(value, (datetime, date)):
            if self.dt_format == 'iso8601':
                return value.isoformat()
            else:
                return value.strftime(self.dt_format)
        
        # For any other type, convert to string
        return str(value)
# Flask-RESTful response fields
molecule_fields = {
    'id': fields.String,
    'cid': fields.Integer,
    'name': fields.String,
    'molecular_formula': fields.String,
    'smiles': fields.String,
    'pubchem_link': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601'),
    'updated_at': FlexibleDateTime(dt_format='iso8601'),
    'properties': fields.Raw
}

mixture_fields = {
    'id': fields.String,
    'name': fields.String,
    'description': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601'),
    'updated_at': FlexibleDateTime(dt_format='iso8601'),
    'created_by': fields.String,
    'components': fields.Raw
}

prediction_fields = {
    'id': fields.String,
    'mixture_id': fields.String,
    'property_type_id': fields.String,
    'property_name': fields.String,
    'numeric_value': fields.Float,
    'text_value': fields.String,
    'boolean_value': fields.Boolean,
    'confidence': fields.Float,
    'calculation_method': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601')
}

experiment_fields = {
    'id': fields.String,
    'mixture_id': fields.String,
    'property_type_id': fields.String,
    'property_name': fields.String,
    'numeric_value': fields.Float,
    'text_value': fields.String,
    'boolean_value': fields.Boolean,
    'experimental_conditions': fields.String,
    'date_performed': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601')
}

comparison_fields = {
    'prediction': fields.Raw,
    'experiment': fields.Raw,
    'difference': fields.Float,
    'percent_error': fields.Float
}

# RDKit-related fields
molecular_property_fields = {
    'hydrogen_bonding': fields.Raw,
    'logp': fields.Float,
    'tpsa': fields.Float,
    'molecular_properties': fields.Raw,
    'functional_groups': fields.Raw,
    'permeability': fields.Raw,
    'smiles': fields.String,
    'inchi': fields.String,
    'inchi_key': fields.String
}

molecule_visualization_fields = {
    'svg': fields.String,
    'width': fields.Integer,
    'height': fields.Integer
}

substructure_search_fields = {
    'match': fields.Boolean,
    'match_count': fields.Integer,
    'matches': fields.Raw,
    'visualization': fields.String
}

similarity_fields = {
    'tanimoto': fields.Float,
    'dice': fields.Float,
    'fingerprint_type': fields.String
}

# Toxicity-related fields
toxicity_data_fields = {
    'id': fields.String,
    'molecule_id': fields.String,
    'assay_id': fields.String,
    'activity_value': fields.Float,
    'activity_unit': fields.String,
    'activity_type': fields.String,
    'hit_call': fields.Boolean,
    'significance': fields.Float,
    'reliability_score': fields.Float,
    'data_quality_comment': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601')
}

toxicity_assay_fields = {
    'id': fields.String,
    'source_id': fields.String,
    'assay_name': fields.String,
    'assay_id': fields.String,
    'description': fields.String,
    'assay_type': fields.String,
    'organism': fields.String,
    'tissue': fields.String,
    'cell_line': fields.String,
    'assay_target': fields.String,
    'toxicological_endpoint': fields.String
}

toxicity_endpoint_fields = {
    'id': fields.String,
    'name': fields.String,
    'description': fields.String,
    'category': fields.String
}

toxicity_score_fields = {
    'id': fields.String,
    'molecule_id': fields.String,
    'score_type': fields.String,
    'score_value': fields.Float,
    'confidence': fields.Float,
    'method_id': fields.String,
    'endpoint_id': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601')
}

unified_score_fields = {
    'score': fields.Float,
    'application_context': fields.String,
    'component_scores': fields.Raw,
    'confidence': fields.Float,
    'details': fields.Raw
}

# Marshmallow validation schemas
class MoleculeSchema(Schema):
    """
    Schema for creating or updating a molecule.
    
    This schema validates the data for a molecule, including its basic
    information and structural representation. It ensures that all required
    fields are present and that values meet validation requirements.
    
    Attributes:
        cid: PubChem Compound ID (optional for user-created molecules)
        name: Name of the molecule (required)
        molecular_formula: Chemical formula of the molecule (required)
        smiles: SMILES notation representing molecular structure (required)
        pubchem_link: Link to PubChem page (generated automatically for imported molecules)
    
    Validation:
        - Name must be between 1 and 200 characters
        - SMILES notation must be a valid string
        - Formula must follow chemical formula conventions
    """
    id = ma_fields.UUID(
        description="UUID of the molecule in the database",
        example="123e4567-e89b-12d3-a456-426614174000",
        dump_only=True  # Read-only field, not required for creation/updates
    )
    cid = ma_fields.Integer(
        description="PubChem Compound ID",
        example=2244,
        allow_none=True
    )
    name = ma_fields.String(
        required=True,
        description="Name of the molecule",
        example="Aspirin",
        validate=validate.Length(min=1, max=200)
    )
    molecular_formula = ma_fields.String(
        required=True,
        description="Chemical formula of the molecule",
        example="C9H8O4"
    )
    smiles = ma_fields.String(
        required=True,
        description="SMILES notation representing molecular structure",
        example="CC(=O)OC1=CC=CC=C1C(=O)O"
    )
    pubchem_link = ma_fields.String(
        description="Link to PubChem page",
        example="https://pubchem.ncbi.nlm.nih.gov/compound/2244",
        dump_only=True  # Read-only field, generated automatically
    )
    created_at = ma_fields.DateTime(
        description="Timestamp when the molecule was created",
        dump_only=True  # Read-only field
    )
    updated_at = ma_fields.DateTime(
        description="Timestamp when the molecule was last updated",
        dump_only=True  # Read-only field
    )
    properties = ma_fields.Raw(
        description="Molecular properties",
        dump_only=True  # Read-only field
    )

    @post_load
    def validate_molecule(self, data, **kwargs):
        """
        Perform additional validation on molecule data.
        
        Args:
            data: The validated data
            **kwargs: Additional arguments passed to the method
            
        Returns:
            dict: The validated data
            
        Raises:
            ValidationError: If validation fails
        """
        # Additional validation could be added here, such as:
        # - Validating SMILES format
        # - Checking formula consistency
        # - Validating against chemical databases
        return data

class MoleculeComponentSchema(Schema):
    """
    Schema for a molecule component in a mixture.
    
    This schema validates the data for a single component in a mixture,
    including the molecule reference, concentration, and unit. It ensures
    that all required fields are present and that values meet validation
    requirements.
    
    Attributes:
        id: UUID of the component in the database (read-only)
        mixture_id: UUID of the mixture this component belongs to (read-only)
        molecule_id: UUID of the molecule in the database (required)
        concentration: Concentration value (must be positive) (required)
        concentration_unit: Unit of concentration (required)
        created_by: UUID of the user who created the component (read-only)
        created_at: Timestamp when the component was created (read-only)
    
    Validation:
        - Concentration must be positive
        - Concentration unit must be one of the allowed values
    
    Response Fields:
        When components are returned as part of a mixture response, they may include
        additional molecule details such as name, formula, and properties.
    """
    id = ma_fields.UUID(
        description="UUID of the component in the database",
        example="123e4567-e89b-12d3-a456-426614174000",
        dump_only=True  # Read-only field
    )
    mixture_id = ma_fields.UUID(
        description="UUID of the mixture this component belongs to",
        example="123e4567-e89b-12d3-a456-426614174000",
        dump_only=True  # Read-only field
    )
    molecule_id = ma_fields.UUID(
        required=True,
        description="UUID of the molecule in the database",
        example="123e4567-e89b-12d3-a456-426614174000"
    )
    concentration = ma_fields.Float(
        required=True,
        validate=validate.Range(min=0),
        description="Concentration value (must be positive)",
        example=10.5
    )
    concentration_unit = ma_fields.String(
        required=True,
        description="Unit of concentration (e.g., 'mg/mL', '%w/v')",
        example="mg/mL",
        validate=validate.OneOf(["mg/mL", "%w/v", "%v/v", "mM", "µM"])
    )
    molecule = ma_fields.Nested("MoleculeSchema", dump_only=True)
    created_by = ma_fields.UUID(
        description="UUID of the user who created the component",
        dump_only=True  # Read-only field
    )
    created_at = ma_fields.DateTime(
        description="Timestamp when the component was created",
        dump_only=True  # Read-only field
    )

class MixtureSchema(Schema):
    """
    Schema for creating or updating a mixture.
    
    This schema validates the data for a mixture, including its basic
    information and a list of components. It ensures that all required
    fields are present and that values meet validation requirements.
    
    Attributes:
        id: UUID of the mixture (read-only)
        name: Name of the mixture (required)
        description: Description of the mixture (optional)
        components: List of components in the mixture (required)
        created_at: Timestamp when the mixture was created (read-only)
        updated_at: Timestamp when the mixture was last updated (read-only)
        created_by: UUID of the user who created the mixture (read-only)
    
    Validation:
        - Name must be between 1 and 100 characters
        - Description must be at most 500 characters
        - Mixture must have at least one component
    
    Response Fields:
        The response includes all the fields above, with components expanded to include
        molecule details such as name, formula, and properties.
    """
    id = ma_fields.UUID(
        description="UUID of the mixture in the database",
        example="123e4567-e89b-12d3-a456-426614174000",
        dump_only=True  # Read-only field, not required for creation/updates
    )
    name = ma_fields.String(
        required=True,
        description="Name of the mixture",
        example="Cryoprotective Solution A",
        validate=validate.Length(min=1, max=100)
    )
    description = ma_fields.String(
        description="Description of the mixture",
        example="Standard cryoprotective solution for cell preservation",
        validate=validate.Length(max=500)
    )
    components = ma_fields.List(
        ma_fields.Nested("MoleculeComponentSchema"),
        required=True,
        description="List of components in the mixture",
        validate=validate.Length(min=1)
    )
    created_at = ma_fields.DateTime(
        description="Timestamp when the mixture was created",
        dump_only=True  # Read-only field
    )
    updated_at = ma_fields.DateTime(
        description="Timestamp when the mixture was last updated",
        dump_only=True  # Read-only field
    )
    created_by = ma_fields.UUID(
        description="UUID of the user who created the mixture",
        dump_only=True  # Read-only field
    )

    @post_load
    def validate_components(self, data, **kwargs):
        """
        Validate that the mixture has at least one component.
        
        Args:
            data: The validated data
            **kwargs: Additional arguments passed to the method
            
        Returns:
            dict: The validated data
            
        Raises:
            ValidationError: If the mixture has no components
        """
        if not data.get('components'):
            raise ValidationError('Mixture must have at least one component')
        
        # Additional validation could be added here, such as:
        # - Validating that all molecule_ids exist in the database
        # - Checking for duplicate molecules in the components list
        # - Validating concentration ranges for specific concentration units
        
        return data

class PropertyValueSchema(Schema):
    """
    Schema for a property value.
    
    This schema validates the data for a property value, including the property name,
    value, and data type. It ensures that the value matches the specified data type.
    """
    property_name = ma_fields.String(
        required=True,
        description="Name of the property",
        example="viscosity"
    )
    value = ma_fields.Raw(
        required=True,
        description="Value of the property (can be numeric, text, or boolean)",
        example=1.25
    )
    data_type = ma_fields.String(
        validate=validate.OneOf(['numeric', 'text', 'boolean']),
        description="Data type of the value",
        example="numeric"
    )

    @post_load
    def validate_value_type(self, data, **kwargs):
        """Validate that the value matches the specified data type."""
        value = data.get('value')
        data_type = data.get('data_type')
        
        if data_type == 'numeric' and not isinstance(value, (int, float)):
            raise ValidationError('Value must be a number for numeric data type')
        elif data_type == 'text' and not isinstance(value, str):
            raise ValidationError('Value must be a string for text data type')
        elif data_type == 'boolean' and not isinstance(value, bool):
            raise ValidationError('Value must be a boolean for boolean data type')
        
        return data

class PredictionSchema(Schema):
    """
    Schema for creating a prediction.
    
    This schema validates the data for a prediction, including the property name,
    predicted value, confidence level, and calculation method used. It ensures
    that all required fields are present and that values meet validation
    requirements.
    
    Attributes:
        property_name: Name of the property being predicted (required)
        value: Predicted value (can be numeric, text, or boolean) (required)
        confidence: Confidence level of the prediction (0-1) (required)
        calculation_method: Method used to calculate the prediction (required)
    
    Validation:
        - Confidence must be between 0 and 1
    """
    property_name = ma_fields.String(
        required=True,
        description="Name of the property being predicted",
        example="freezing_point"
    )
    value = ma_fields.Raw(
        required=True,
        description="Predicted value (can be numeric, text, or boolean)",
        example=-20.5
    )
    confidence = ma_fields.Float(
        required=True,
        validate=validate.Range(min=0, max=1),
        description="Confidence level of the prediction (0-1)",
        example=0.85
    )
    calculation_method = ma_fields.String(
        required=True,
        description="Method used to calculate the prediction",
        example="molecular_dynamics"
    )

class ExperimentSchema(Schema):
    """
    Schema for recording an experiment.
    
    This schema validates the data for an experiment, including the property name,
    measured value, experimental conditions, and the date the experiment was performed.
    """
    property_name = ma_fields.String(
        required=True,
        description="Name of the property being measured",
        example="glass_transition_temperature"
    )
    value = ma_fields.Raw(
        required=True,
        description="Measured value (can be numeric, text, or boolean)",
        example=-45.2
    )
    experimental_conditions = ma_fields.String(
        description="Description of the experimental conditions",
        example="Measured using differential scanning calorimetry at a rate of 10°C/min"
    )
    date_performed = ma_fields.Date(
        required=True,
        description="Date when the experiment was performed (YYYY-MM-DD)",
        example="2024-04-15"
    )

class ComparisonQuerySchema(Schema):
    """
    Schema for querying a comparison.
    
    This schema validates the query parameters for comparing predictions with
    experimental results for a specific property. It ensures that all required
    fields are present for the comparison operation.
    
    Attributes:
        property_name: Name of the property to compare (required)
    
    Usage:
        This schema is used when requesting a comparison between predicted and
        experimental values for a specific property of a mixture.
    """
    property_name = ma_fields.String(
        required=True,
        description="Name of the property to compare",
        example="thermal_conductivity"
    )

# Note: RDKit-related schemas have been moved to api/rdkit_schemas.py

# Base Model class for common functionality
class BaseModel(Generic[T]):
    """
    Base model class with common database operations.
    
    This class provides a foundation for all database models, implementing
    common operations like create, read, update, and delete. It uses the
    Supabase client for database interactions and includes error handling.
    
    Attributes:
        table_name (str): Name of the database table (must be set by subclasses)
        id_field (str): Name of the ID field (default: 'id')
    """
    
    table_name: str = None
    id_field: str = 'id'
    
    @classmethod
    def get_supabase(cls):
        """
        Get the Supabase client.
        
        Returns:
            SupabaseClient: The configured Supabase client
        """
        return get_supabase_client()
    
    @classmethod
    def get_table(cls):
        """
        Get the Supabase table reference.
        
        Returns:
            TableReference: The Supabase table reference
            
        Raises:
            ValueError: If table_name is not defined in the subclass
        """
        if not cls.table_name:
            raise ValueError(f"table_name not defined for {cls.__name__}")
        return cls.get_supabase().table(cls.table_name)
    
    @classmethod
    def create(cls, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create a new record in the database.

        Args:
            data: Dictionary containing the data to insert

        Returns:
            Dict[str, Any]: The created record
            
        Raises:
            Exception: If an error occurs during the database operation
            
        Performance Optimization:
            - Automatically adds user ID for tracking
            - Uses standardized error handling for consistent responses
        """
        # Add created_by if not present
        if 'created_by' not in data:
            user_id = get_user_id()
            if user_id:
                data['created_by'] = user_id

        try:
            response = cls.get_table().insert(data).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Creating {cls.__name__}"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error creating {cls.__name__}: {error_message}")
            return response.data[0] if response.data else None
        except Exception as e:
            context = f"Creating {cls.__name__}"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error creating {cls.__name__}: {str(e)}")
    
    @classmethod
    def get(cls, id_value: str) -> Dict[str, Any]:
        """
        Get a record by ID.

        Args:
            id_value: The ID of the record to retrieve

        Returns:
            Dict[str, Any]: The retrieved record or None if not found
            
        Raises:
            Exception: If an error occurs during the database operation
            
        Performance Optimization:
            - Uses standardized error handling for consistent responses
        """
        try:
            response = cls.get_table().select('*').eq(cls.id_field, id_value).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Getting {cls.__name__} with ID {id_value}"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error getting {cls.__name__}: {error_message}")
            return response.data[0] if response.data else None
        except Exception as e:
            context = f"Getting {cls.__name__} with ID {id_value}"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error getting {cls.__name__}: {str(e)}")
    
    @classmethod
    def get_all(cls, limit: int = 100, offset: int = 0) -> List[Dict[str, Any]]:
        """
        Get all records with pagination.

        Args:
            limit: Maximum number of records to return
            offset: Number of records to skip

        Returns:
            List of records
        """
        try:
            response = cls.get_table().select('*').range(offset, offset + limit - 1).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Getting all {cls.__name__} records"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error getting all {cls.__name__}: {error_message}")
            return response.data
        except Exception as e:
            context = f"Getting all {cls.__name__} records"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error getting all {cls.__name__}: {str(e)}")
    
    @classmethod
    def update(cls, id_value: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update a record by ID.

        Args:
            id_value: The ID of the record to update
            data: Dictionary containing the data to update

        Returns:
            The updated record
        """
        try:
            response = cls.get_table().update(data).eq(cls.id_field, id_value).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Updating {cls.__name__} with ID {id_value}"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error updating {cls.__name__}: {error_message}")
            return response.data[0] if response.data else None
        except Exception as e:
            context = f"Updating {cls.__name__} with ID {id_value}"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error updating {cls.__name__}: {str(e)}")
    
    @classmethod
    def delete(cls, id_value: str) -> bool:
        """
        Delete a record by ID.

        Args:
            id_value: The ID of the record to delete

        Returns:
            True if successful, False otherwise
        """
        try:
            response = cls.get_table().delete().eq(cls.id_field, id_value).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Deleting {cls.__name__} with ID {id_value}"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error deleting {cls.__name__}: {error_message}")
            return True
        except Exception as e:
            context = f"Deleting {cls.__name__} with ID {id_value}"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error deleting {cls.__name__}: {str(e)}")
    
    @classmethod
    def filter(cls, filters: Dict[str, Any], limit: int = 100, offset: int = 0) -> List[Dict[str, Any]]:
        """
        Filter records by multiple criteria.

        Args:
            filters: Dictionary of field-value pairs to filter by
            limit: Maximum number of records to return
            offset: Number of records to skip

        Returns:
            List of matching records
        """
        try:
            query = cls.get_table().select('*')
            for field, value in filters.items():
                query = query.eq(field, value)
            query = query.range(offset, offset + limit - 1)
            response = query.execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Filtering {cls.__name__} records"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error filtering {cls.__name__}: {error_message}")
            return response.data
        except Exception as e:
            context = f"Filtering {cls.__name__} records"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error filtering {cls.__name__}: {str(e)}")
    
    @classmethod
    def count(cls, filters: Dict[str, Any] = None) -> int:
        """
        Count records, optionally filtered.

        Args:
            filters: Optional dictionary of field-value pairs to filter by

        Returns:
            Count of matching records
        """
        try:
            query = cls.get_table().select('id', count='exact')
            if filters:
                for field, value in filters.items():
                    query = query.eq(field, value)
            response = query.execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Counting {cls.__name__} records"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error counting {cls.__name__}: {error_message}")
            return response.count if hasattr(response, 'count') else len(response.data)
        except Exception as e:
            context = f"Counting {cls.__name__} records"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error counting {cls.__name__}: {str(e)}")
    
    @classmethod
    def exists(cls, id_value: str) -> bool:
        """
        Check if a record exists by ID.

        Args:
            id_value: The ID to check

        Returns:
            True if the record exists, False otherwise
        """
        try:
            response = cls.get_table().select('id').eq(cls.id_field, id_value).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Checking existence of {cls.__name__} with ID {id_value}"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error checking existence of {cls.__name__}: {error_message}")
            return len(response.data) > 0
        except Exception as e:
            context = f"Checking existence of {cls.__name__} with ID {id_value}"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error checking existence of {cls.__name__}: {str(e)}")


# Molecule model
class Molecule(BaseModel):
    """Model for molecules in the database."""
    
    table_name = 'molecules'
    
    @classmethod
    def create_from_pubchem(cls, cid: int) -> Dict[str, Any]:
        """
        Import a molecule from PubChem.

        Args:
            cid: PubChem Compound ID

        Returns:
            The imported molecule record
        """
        user_id = get_user_id()
        try:
            response = cls.get_supabase().rpc(
                "import_molecule_from_pubchem",
                {"p_cid": cid, "p_user_id": user_id}
            ).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                raise Exception(f"Error importing molecule from PubChem: {error_message}")
            molecule_id = response.data
            return cls.get(molecule_id)
        except Exception as e:
            raise Exception(f"Error importing molecule from PubChem: {str(e)}")
    
    @classmethod
    def get_with_properties(cls, id_value: str) -> Dict[str, Any]:
        """
        Get a molecule with all its properties.

        Args:
            id_value: The ID of the molecule

        Returns:
            Molecule with properties
        """
        try:
            response = cls.get_supabase().table("molecule_with_properties").select("*").eq("id", id_value).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                raise Exception(f"Error getting molecule with properties: {error_message}")
            return response.data[0] if response.data else None
        except Exception as e:
            raise Exception(f"Error getting molecule with properties: {str(e)}")
    
    @classmethod
    def search_by_name(cls, name: str, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search molecules by name.

        Args:
            name: Name to search for (case-insensitive, partial match)
            limit: Maximum number of results to return

        Returns:
            List of matching molecules
        """
        try:
            response = cls.get_table().select('*').ilike('name', f'%{name}%').limit(limit).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                raise Exception(f"Error searching molecules by name: {error_message}")
            return response.data
        except Exception as e:
            raise Exception(f"Error searching molecules by name: {str(e)}")
    
    @classmethod
    def search_by_formula(cls, formula: str, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search molecules by molecular formula.

        Args:
            formula: Formula to search for (case-sensitive, exact match)
            limit: Maximum number of results to return

        Returns:
            List of matching molecules
        """
        try:
            response = cls.get_table().select('*').eq('molecular_formula', formula).limit(limit).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                raise Exception(f"Error searching molecules by formula: {error_message}")
            return response.data
        except Exception as e:
            raise Exception(f"Error searching molecules by formula: {str(e)}")
    
    @classmethod
    def get_by_cid(cls, cid: int) -> Dict[str, Any]:
        """
        Get a molecule by PubChem CID.

        Args:
            cid: PubChem Compound ID

        Returns:
            Molecule record or None if not found
        """
        try:
            response = cls.get_table().select('*').eq('cid', cid).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                raise Exception(f"Error getting molecule by CID: {error_message}")
            return response.data[0] if response.data else None
        except Exception as e:
            raise Exception(f"Error getting molecule by CID: {str(e)}")


# Property Type model
class PropertyType(BaseModel):
    """Model for property types in the database."""
    
    table_name = 'property_types'
    
    @classmethod
    def get_by_name(cls, name: str) -> Dict[str, Any]:
        """
        Get a property type by name.
        
        Args:
            name: Name of the property type
            
        Returns:
            Property type record or None if not found
        """
        response = cls.get_table().select('*').eq('name', name).execute()
        
        if response.error:
            raise Exception(f"Error getting property type by name: {response.error}")
        
        return response.data[0] if response.data else None


# Molecular Property model
class MolecularProperty(BaseModel):
    """Model for molecular properties in the database."""
    
    table_name = 'molecular_properties'
    
    @classmethod
    def add_property(cls, molecule_id: str, property_name: str, value: Any) -> Dict[str, Any]:
        """
        Add a property to a molecule.
        
        Args:
            molecule_id: ID of the molecule
            property_name: Name of the property
            value: Value of the property
            
        Returns:
            The created property record
        """
        # Get property type
        property_type = PropertyType.get_by_name(property_name)
        if not property_type:
            raise ValueError(f"Property type '{property_name}' not found")
        
        # Prepare property data
        property_data = {
            "molecule_id": molecule_id,
            "property_type_id": property_type["id"],
            "created_by": get_user_id()
        }
        
        # Set the appropriate value field based on data type
        if property_type["data_type"] == "numeric":
            property_data["numeric_value"] = float(value) if value is not None else None
        elif property_type["data_type"] == "text":
            property_data["text_value"] = str(value) if value is not None else None
        elif property_type["data_type"] == "boolean":
            property_data["boolean_value"] = bool(value) if value is not None else None
        else:
            raise ValueError(f"Unknown data type '{property_type['data_type']}'")
        
        # Insert or update property
        response = cls.get_table().upsert(
            property_data,
            on_conflict="molecule_id,property_type_id"
        ).execute()
        
        if response.error:
            raise Exception(f"Error adding property: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def add_multiple_properties(cls, molecule_id: str, properties: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Add multiple properties to a molecule.
        
        Args:
            molecule_id: ID of the molecule
            properties: Dictionary of property name-value pairs
            
        Returns:
            List of created property records
        """
        # Get property types
        response = cls.get_supabase().table("property_types").select("id, name, data_type").execute()
        if response.error:
            raise Exception(f"Error fetching property types: {response.error}")
        
        property_types = response.data
        
        # Get the current user ID
        user_id = get_user_id()
        
        # Prepare property inserts
        property_inserts = []
        
        for property_name, value in properties.items():
            property_type = next((pt for pt in property_types if pt["name"] == property_name), None)
            if not property_type:
                continue
            
            property_insert = {
                "molecule_id": molecule_id,
                "property_type_id": property_type["id"],
                "created_by": user_id
            }
            
            # Set the appropriate value field based on data type
            if property_type["data_type"] == "numeric":
                property_insert["numeric_value"] = float(value) if value is not None else None
            elif property_type["data_type"] == "text":
                property_insert["text_value"] = str(value) if value is not None else None
            elif property_type["data_type"] == "boolean":
                property_insert["boolean_value"] = bool(value) if value is not None else None
            else:
                continue
            
            property_inserts.append(property_insert)
        
        # Insert properties
        if not property_inserts:
            return []
            
        response = cls.get_table().upsert(
            property_inserts,
            on_conflict="molecule_id,property_type_id"
        ).execute()
        
        if response.error:
            raise Exception(f"Error adding properties: {response.error}")
        
        return response.data
    
    @classmethod
    def get_property(cls, molecule_id: str, property_name: str) -> Dict[str, Any]:
        """
        Get a specific property of a molecule.
        
        Args:
            molecule_id: ID of the molecule
            property_name: Name of the property
            
        Returns:
            Property record or None if not found
        """
        # Get property type
        property_type = PropertyType.get_by_name(property_name)
        if not property_type:
            raise ValueError(f"Property type '{property_name}' not found")
        
        response = cls.get_table().select('*').eq('molecule_id', molecule_id).eq('property_type_id', property_type["id"]).execute()
        
        if response.error:
            raise Exception(f"Error getting property: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def get_properties(cls, molecule_id: str) -> List[Dict[str, Any]]:
        """
        Get all properties of a molecule.
        
        Args:
            molecule_id: ID of the molecule
            
        Returns:
            List of property records
        """
        response = cls.get_table().select('*').eq('molecule_id', molecule_id).execute()
        
        if response.error:
            raise Exception(f"Error getting properties: {response.error}")
        
        return response.data


# Mixture model
class Mixture(BaseModel):
    """Model for mixtures in the database."""
    
    table_name = 'mixtures'
    
    @classmethod
    def create_with_components(cls, name: str, description: str, components: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Create a mixture with components.
        
        Args:
            name: Name of the mixture
            description: Description of the mixture
            components: List of dicts with molecule_id, concentration, concentration_unit
            
        Returns:
            The created mixture record
        """
        # Get the current user ID
        user_id = get_user_id()
        
        # Create mixture
        mixture_data = {
            "name": name,
            "description": description,
            "created_by": user_id
        }
        
        response = cls.get_table().insert(mixture_data).execute()
        
        if response.error:
            raise Exception(f"Error creating mixture: {response.error}")
        
        mixture_id = response.data[0]["id"]
        
        # Add components
        component_inserts = [{
            "mixture_id": mixture_id,
            "molecule_id": comp["molecule_id"],
            "concentration": comp["concentration"],
            "concentration_unit": comp["concentration_unit"],
            "created_by": user_id
        } for comp in components]
        
        response = cls.get_supabase().table("mixture_components").insert(component_inserts).execute()
        
        if response.error:
            raise Exception(f"Error adding mixture components: {response.error}")
        
        return cls.get_with_components(mixture_id)
    
    @classmethod
    def get_with_components(cls, id_value: str) -> Dict[str, Any]:
        """
        Get a mixture with all its components.
        
        Args:
            id_value: The ID of the mixture
            
        Returns:
            Mixture with components
        """
        response = cls.get_supabase().table("mixture_with_components").select("*").eq("id", id_value).execute()
        
        if response.error:
            raise Exception(f"Error getting mixture with components: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def update_with_components(cls, id_value: str, data: Dict[str, Any], components: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Update a mixture and its components.
        
        Args:
            id_value: The ID of the mixture to update
            data: Dictionary containing the mixture data to update
            components: List of dicts with molecule_id, concentration, concentration_unit
            
        Returns:
            The updated mixture record
        """
        # Update mixture
        response = cls.get_table().update(data).eq(cls.id_field, id_value).execute()
        
        if response.error:
            raise Exception(f"Error updating mixture: {response.error}")
        
        # Delete existing components
        response = cls.get_supabase().table("mixture_components").delete().eq("mixture_id", id_value).execute()
        
        if response.error:
            raise Exception(f"Error deleting mixture components: {response.error}")
        
        # Add new components
        user_id = get_user_id()
        component_inserts = [{
            "mixture_id": id_value,
            "molecule_id": comp["molecule_id"],
            "concentration": comp["concentration"],
            "concentration_unit": comp["concentration_unit"],
            "created_by": user_id
        } for comp in components]
        
        response = cls.get_supabase().table("mixture_components").insert(component_inserts).execute()
        
        if response.error:
            raise Exception(f"Error adding mixture components: {response.error}")
        
        return cls.get_with_components(id_value)
    
    @classmethod
    def search_by_name(cls, name: str, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search mixtures by name.
        
        Args:
            name: Name to search for (case-insensitive, partial match)
            limit: Maximum number of results to return
            
        Returns:
            List of matching mixtures
        """
        response = cls.get_table().select('*').ilike('name', f'%{name}%').limit(limit).execute()
        
        if response.error:
            raise Exception(f"Error searching mixtures by name: {response.error}")
        
        return response.data
    
    @classmethod
    def get_mixtures_with_molecule(cls, molecule_id: str) -> List[Dict[str, Any]]:
        """
        Get all mixtures containing a specific molecule.
        
        Args:
            molecule_id: ID of the molecule
            
        Returns:
            List of mixtures
        """
        # First get mixture IDs from mixture_components
        response = cls.get_supabase().table("mixture_components").select("mixture_id").eq("molecule_id", molecule_id).execute()
        
        if response.error:
            raise Exception(f"Error getting mixtures with molecule: {response.error}")
        
        if not response.data:
            return []
        
        # Then get the mixtures
        mixture_ids = [item["mixture_id"] for item in response.data]
        response = cls.get_table().select('*').in_('id', mixture_ids).execute()
        
        if response.error:
            raise Exception(f"Error getting mixtures: {response.error}")
        
        return response.data
    
    @classmethod
    def calculate_score(cls, id_value: str) -> float:
        """
        Calculate the score for a mixture.
        
        Args:
            id_value: The ID of the mixture
            
        Returns:
            The calculated score
        """
        response = cls.get_supabase().rpc(
            "calculate_mixture_score",
            {"p_mixture_id": id_value}
        ).execute()
        
        if response.error:
            raise Exception(f"Error calculating mixture score: {response.error}")
        
        return response.data


# Calculation Method model
class CalculationMethod(BaseModel):
    """Model for calculation methods in the database."""
    
    table_name = 'calculation_methods'
    
    @classmethod
    def get_by_name(cls, name: str) -> Dict[str, Any]:
        """
        Get a calculation method by name.
        
        Args:
            name: Name of the calculation method
            
        Returns:
            Calculation method record or None if not found
        """
        response = cls.get_table().select('*').eq('name', name).execute()
        
        if response.error:
            raise Exception(f"Error getting calculation method by name: {response.error}")
        
        return response.data[0] if response.data else None


# Prediction model
class Prediction(BaseModel):
    """Model for predictions in the database."""
    
    table_name = 'predictions'
    
    @classmethod
    def add_prediction(cls, mixture_id: str, property_name: str, value: Any, confidence: float, method_name: str) -> Dict[str, Any]:
        """
        Add a prediction for a mixture.
        
        Args:
            mixture_id: ID of the mixture
            property_name: Name of the property being predicted
            value: Predicted value
            confidence: Confidence level (0-1)
            method_name: Name of the calculation method
        
        Returns:
            The inserted prediction record
        """
        # Get property type ID
        property_type = PropertyType.get_by_name(property_name)
        if not property_type:
            raise ValueError(f"Property type '{property_name}' not found")
        
        # Get calculation method ID
        calculation_method = CalculationMethod.get_by_name(method_name)
        if not calculation_method:
            raise ValueError(f"Calculation method '{method_name}' not found")
        
        # Get the current user ID
        user_id = get_user_id()
        
        # Prepare prediction insert
        prediction = {
            "mixture_id": mixture_id,
            "property_type_id": property_type["id"],
            "calculation_method_id": calculation_method["id"],
            "confidence": confidence,
            "created_by": user_id
        }
        
        # Set the appropriate value field based on data type
        if property_type["data_type"] == "numeric":
            prediction["numeric_value"] = float(value) if value is not None else None
        elif property_type["data_type"] == "text":
            prediction["text_value"] = str(value) if value is not None else None
        elif property_type["data_type"] == "boolean":
            prediction["boolean_value"] = bool(value) if value is not None else None
        else:
            raise ValueError(f"Unknown data type '{property_type['data_type']}'")
        
        # Insert prediction
        response = cls.get_table().upsert(
            prediction,
            on_conflict="mixture_id,property_type_id,calculation_method_id"
        ).execute()
        
        if response.error:
            raise Exception(f"Error adding prediction: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def get_predictions_for_mixture(cls, mixture_id: str) -> List[Dict[str, Any]]:
        """
        Get all predictions for a mixture.
        
        Args:
            mixture_id: ID of the mixture
            
        Returns:
            List of prediction records
        """
        response = cls.get_table().select('*').eq('mixture_id', mixture_id).execute()
        
        if response.error:
            raise Exception(f"Error getting predictions for mixture: {response.error}")
        
        return response.data
    
    @classmethod
    def get_prediction(cls, mixture_id: str, property_name: str, method_name: str) -> Dict[str, Any]:
        """
        Get a specific prediction for a mixture.
        
        Args:
            mixture_id: ID of the mixture
            property_name: Name of the property
            method_name: Name of the calculation method
            
        Returns:
            Prediction record or None if not found
        """
        # Get property type
        property_type = PropertyType.get_by_name(property_name)
        if not property_type:
            raise ValueError(f"Property type '{property_name}' not found")
        
        # Get calculation method
        calculation_method = CalculationMethod.get_by_name(method_name)
        if not calculation_method:
            raise ValueError(f"Calculation method '{method_name}' not found")
        
        response = cls.get_table().select('*').eq('mixture_id', mixture_id).eq('property_type_id', property_type["id"]).eq('calculation_method_id', calculation_method["id"]).execute()
        
        if response.error:
            raise Exception(f"Error getting prediction: {response.error}")
        
        return response.data[0] if response.data else None


# Experiment model
class Experiment(BaseModel):
    """Model for experiments in the database."""
    
    table_name = 'experiments'
    
    @classmethod
    def record_experiment(cls, mixture_id: str, property_name: str, value: Any, conditions: str, experiment_date: Union[str, date]) -> Dict[str, Any]:
        """
        Record experimental results.
        
        Args:
            mixture_id: ID of the mixture
            property_name: Name of the property being measured
            value: Measured value
            conditions: Experimental conditions
            experiment_date: Date of experiment (YYYY-MM-DD or date object)
            
        Returns:
            The inserted experiment record
        """
        # Get property type ID
        property_type = PropertyType.get_by_name(property_name)
        if not property_type:
            raise ValueError(f"Property type '{property_name}' not found")
        
        # Get the current user ID
        user_id = get_user_id()
        
        # Convert date to string if it's a date object
        if isinstance(experiment_date, date):
            experiment_date = experiment_date.isoformat()
        
        # Prepare experiment insert
        experiment = {
            "mixture_id": mixture_id,
            "property_type_id": property_type["id"],
            "experimental_conditions": conditions,
            "date_performed": experiment_date,
            "created_by": user_id
        }
        
        # Set the appropriate value field based on data type
        if property_type["data_type"] == "numeric":
            experiment["numeric_value"] = float(value) if value is not None else None
        elif property_type["data_type"] == "text":
            experiment["text_value"] = str(value) if value is not None else None
        elif property_type["data_type"] == "boolean":
            experiment["boolean_value"] = bool(value) if value is not None else None
        else:
            raise ValueError(f"Unknown data type '{property_type['data_type']}'")
        
        # Insert experiment
        response = cls.get_table().insert(experiment).execute()
        
        if response.error:
            raise Exception(f"Error recording experiment: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def get_experiments_for_mixture(cls, mixture_id: str) -> List[Dict[str, Any]]:
        """
        Get all experiments for a mixture.
        
        Args:
            mixture_id: ID of the mixture
            
        Returns:
            List of experiment records
        """
        response = cls.get_table().select('*').eq('mixture_id', mixture_id).execute()
        
        if response.error:
            raise Exception(f"Error getting experiments for mixture: {response.error}")
        
        return response.data
    
    @classmethod
    def get_experiment(cls, mixture_id: str, property_name: str) -> Dict[str, Any]:
        """
        Get the most recent experiment for a mixture and property.
        
        Args:
            mixture_id: ID of the mixture
            property_name: Name of the property
            
        Returns:
            Experiment record or None if not found
        """
        # Get property type
        property_type = PropertyType.get_by_name(property_name)
        if not property_type:
            raise ValueError(f"Property type '{property_name}' not found")
        
        response = cls.get_table().select('*').eq('mixture_id', mixture_id).eq('property_type_id', property_type["id"]).order('date_performed', desc=True).order('created_at', desc=True).limit(1).execute()
        
        if response.error:
            raise Exception(f"Error getting experiment: {response.error}")
        
        return response.data[0] if response.data else None


# Comparison model for comparing predictions with experiments
class Comparison:
    """Utility class for comparing predictions with experiments."""
    
    @staticmethod
    def compare_prediction_with_experiment(mixture_id: str, property_name: str) -> Dict[str, Any]:
        """
        Compare a prediction with an experiment.
        
        Args:
            mixture_id: ID of the mixture
            property_name: Name of the property
            
        Returns:
            Comparison result
        """
        # Get property type
        property_type = PropertyType.get_by_name(property_name)
        if not property_type:
            raise ValueError(f"Property type '{property_name}' not found")
        
        # Call the database function to compare
        response = BaseModel.get_supabase().rpc(
            "compare_prediction_with_experiment",
            {
                "p_mixture_id": mixture_id,
                "p_property_type_id": property_type["id"]
            }
        ).execute()
        
        if response.error:
            raise Exception(f"Error comparing prediction with experiment: {response.error}")
        
        return response.data


# Project-related fields
project_fields = {
    'id': fields.String,
    'name': fields.String,
    'description': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601'),
    'updated_at': FlexibleDateTime(dt_format='iso8601'),
    'created_by': fields.String,
    'is_public': fields.Boolean,
    'experiment_count': fields.Integer
}

project_experiment_fields = {
    'id': fields.String,
    'project_id': fields.String,
    'experiment_id': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601'),
    'notes': fields.String
}

# Project-related schemas
class ProjectSchema(Schema):
    """
    Schema for creating or updating a project.
    
    This schema validates the data for a project, including its name,
    description, and visibility settings.
    """
    name = ma_fields.String(
        required=True,
        description="Name of the project",
        example="Cryoprotectant Optimization Study",
        validate=validate.Length(min=1, max=100)
    )
    description = ma_fields.String(
        description="Description of the project",
        example="Investigation of novel cryoprotectant mixtures for cell preservation",
        validate=validate.Length(max=500)
    )
    is_public = ma_fields.Boolean(
        default=False,
        description="Whether the project is publicly visible",
        example=False
    )

class ToxicityQuerySchema(Schema):
    """
    Schema for toxicity data query parameters.
    
    This schema validates query parameters for toxicity data endpoints,
    ensuring that filters and parameters meet validation requirements.
    
    Attributes:
        assay_id: Optional filter by assay ID
        endpoint: Optional filter by toxicological endpoint
        active_only: Optional filter for active results only
    """
    assay_id = ma_fields.UUID(
        description="Filter by assay ID",
        example="123e4567-e89b-12d3-a456-426614174000",
        required=False
    )
    endpoint = ma_fields.String(
        description="Filter by toxicological endpoint",
        example="developmental_toxicity",
        required=False
    )
    active_only = ma_fields.Boolean(
        description="Return only active results",
        example=True,
        required=False,
        default=False
    )

class UnifiedScoreQuerySchema(Schema):
    """
    Schema for unified score query parameters.
    
    This schema validates query parameters for unified scoring endpoints,
    ensuring that parameters meet validation requirements.
    
    Attributes:
        context: Application context for weighting
        algorithm: Algorithm to use for predictive models
        recalculate: Whether to recalculate scores even if they already exist
    """
    context = ma_fields.String(
        description="Application context for weighting",
        example="general",
        required=False,
        default="general",
        validate=validate.OneOf(["general", "cell_preservation", "organ_preservation",
                                "long_term_storage", "sensitive_tissues"])
    )
    algorithm = ma_fields.String(
        description="Algorithm to use for predictive models",
        example="random_forest",
        required=False,
        default="random_forest",
        validate=validate.OneOf(["random_forest", "neural_network", "ensemble"])
    )
    recalculate = ma_fields.Boolean(
        description="Whether to recalculate scores even if they already exist",
        example=False,
        required=False,
        default=False
    )

class ProjectExperimentSchema(Schema):
    """
    Schema for adding an experiment to a project.
    
    This schema validates the data for associating an experiment with a project,
    including the experiment ID and optional notes.
    """
    experiment_id = ma_fields.UUID(
        required=True,
        description="UUID of the experiment to add to the project",
        example="123e4567-e89b-12d3-a456-426614174000"
    )
    notes = ma_fields.String(
        description="Notes about this experiment in the context of the project",
        example="This experiment confirms the viability of the mixture at -80°C",
        validate=validate.Length(max=1000)
    )


class Project(BaseModel):
    """Model for projects in the database."""
    
    table_name = 'projects'
    
    @classmethod
    def get_with_experiment_count(cls, id_value: str) -> Dict[str, Any]:
        """
        Get a project with the count of experiments.
        
        Args:
            id_value: The ID of the project
            
        Returns:
            Project with experiment count
        """
        response = cls.get_supabase().rpc(
            "get_project_with_experiment_count",
            {"p_project_id": id_value}
        ).execute()
        
        if response.error:
            raise Exception(f"Error getting project with experiment count: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def get_user_projects(cls, user_id: str, limit: int = 100, offset: int = 0) -> List[Dict[str, Any]]:
        """
        Get all projects for a user.
        
        Args:
            user_id: The ID of the user
            limit: Maximum number of projects to return
            offset: Number of projects to skip
            
        Returns:
            List of projects
        """
        response = cls.get_supabase().rpc(
            "get_user_projects",
            {"p_user_id": user_id, "p_limit": limit, "p_offset": offset}
        ).execute()
        
        if response.error:
            raise Exception(f"Error getting user projects: {response.error}")
        
        return response.data
    
    @classmethod
    def get_recent_activity(cls, project_id: str, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Get recent activity for a project.
        
        Args:
            project_id: The ID of the project
            limit: Maximum number of activities to return
            
        Returns:
            List of recent activities
        """
        response = cls.get_supabase().rpc(
            "get_project_activity",
            {"p_project_id": project_id, "p_limit": limit}
        ).execute()
        
        if response.error:
            raise Exception(f"Error getting project activity: {response.error}")
        
        return response.data


class ProjectExperiment(BaseModel):
    """Model for project experiments in the database."""
    
    table_name = 'project_experiments'
    
    @classmethod
    def add_experiment_to_project(cls, project_id: str, experiment_id: str, notes: str = None) -> Dict[str, Any]:
        """
        Add an experiment to a project.
        
        Args:
            project_id: The ID of the project
            experiment_id: The ID of the experiment
            notes: Optional notes about the experiment in this project
            
        Returns:
            The created project experiment record
        """
        user_id = get_user_id()
        
        data = {
            "project_id": project_id,
            "experiment_id": experiment_id,
            "notes": notes,
            "created_by": user_id
        }
        
        response = cls.get_table().insert(data).execute()
        
        if response.error:
            raise Exception(f"Error adding experiment to project: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def get_project_experiments(cls, project_id: str) -> List[Dict[str, Any]]:
        """
        Get all experiments for a project.
        
        Args:
            project_id: The ID of the project
            
        Returns:
            List of experiments
        """
        response = cls.get_supabase().rpc(
            "get_project_experiments",
            {"p_project_id": project_id}
        ).execute()
        
        if response.error:
            raise Exception(f"Error getting project experiments: {response.error}")
        
        return response.data
    
    @classmethod
    def remove_experiment_from_project(cls, project_id: str, experiment_id: str) -> bool:
        """
        Remove an experiment from a project.
        
        Args:
            project_id: The ID of the project
            experiment_id: The ID of the experiment
            
        Returns:
            True if successful, False otherwise
        """
        response = cls.get_table().delete().eq("project_id", project_id).eq("experiment_id", experiment_id).execute()
        
        if response.error:
            raise Exception(f"Error removing experiment from project: {response.error}")
        
        return True

class ToxicityData(BaseModel):
    """Model for toxicity data in the database."""
    
    table_name = "toxicity_data"
    
    @classmethod
    def get_for_molecule(cls, molecule_id: str, assay_id: Optional[str] = None,
                        endpoint: Optional[str] = None, active_only: bool = False) -> List[Dict[str, Any]]:
        """
        Get toxicity data for a specific molecule.
        
        Args:
            molecule_id: UUID of the molecule
            assay_id: Optional filter by assay ID
            endpoint: Optional filter by toxicological endpoint
            active_only: Optional filter for active results only
            
        Returns:
            List of toxicity data records
        """
        try:
            # Start query
            query = cls.get_supabase().table(cls.table_name).select("*").eq("molecule_id", molecule_id)
            
            # Apply filters
            if assay_id:
                query = query.eq("assay_id", assay_id)
                
            if active_only:
                query = query.eq("hit_call", True)
                
            # Execute query
            response = query.execute()
            
            if not response.data:
                return []
                
            # If endpoint filter is provided, we need to join with toxicity_assay
            if endpoint and response.data:
                assay_ids = [item["assay_id"] for item in response.data]
                assay_response = cls.get_supabase().table("toxicity_assay").select("id").in_("id", assay_ids).eq("toxicological_endpoint", endpoint).execute()
                
                if not assay_response.data:
                    return []
                    
                filtered_assay_ids = [item["id"] for item in assay_response.data]
                return [item for item in response.data if item["assay_id"] in filtered_assay_ids]
                
            return response.data
            
        except Exception as e:
            logger.error(f"Error getting toxicity data for molecule {molecule_id}: {str(e)}")
            return []

class ToxicityAssay(BaseModel):
    """Model for toxicity assays in the database."""
    
    table_name = "toxicity_assay"
    
    @classmethod
    def get_all(cls, endpoint: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Get all toxicity assays.
        
        Args:
            endpoint: Optional filter by toxicological endpoint
            
        Returns:
            List of toxicity assays
        """
        try:
            # Start query
            query = cls.get_supabase().table(cls.table_name).select("*")
            
            # Apply filters
            if endpoint:
                query = query.eq("toxicological_endpoint", endpoint)
                
            # Execute query
            response = query.execute()
            
            return response.data or []
            
        except Exception as e:
            logger.error(f"Error getting toxicity assays: {str(e)}")
            return []

class ToxicityEndpoint(BaseModel):
    """Model for toxicity endpoints in the database."""
    
    table_name = "toxicity_endpoint"
    
    @classmethod
    def get_all(cls) -> List[Dict[str, Any]]:
        """
        Get all toxicity endpoints.
        
        Returns:
            List of toxicity endpoints
        """
        try:
            response = cls.get_supabase().table(cls.table_name).select("*").execute()
            return response.data or []
            
        except Exception as e:
            logger.error(f"Error getting toxicity endpoints: {str(e)}")
            return []
            
    @classmethod
    def get_by_ids(cls, endpoint_ids: List[str]) -> List[Dict[str, Any]]:
        """
        Get toxicity endpoints by IDs.
        
        Args:
            endpoint_ids: List of endpoint IDs
            
        Returns:
            List of toxicity endpoints
        """
        try:
            response = cls.get_supabase().table(cls.table_name).select("*").in_("id", endpoint_ids).execute()
            return response.data or []
            
        except Exception as e:
            logger.error(f"Error getting toxicity endpoints by IDs: {str(e)}")
            return []

class ToxicityScore(BaseModel):
    """Model for toxicity scores in the database."""
    
    table_name = "toxicity_score"
    
    @classmethod
    def get_for_molecule(cls, molecule_id: str) -> List[Dict[str, Any]]:
        """
        Get toxicity scores for a specific molecule.
        
        Args:
            molecule_id: UUID of the molecule
            
        Returns:
            List of toxicity score records
        """
        try:
            response = cls.get_supabase().table(cls.table_name).select("*").eq("molecule_id", molecule_id).execute()
            return response.data or []
            
        except Exception as e:
            logger.error(f"Error getting toxicity scores for molecule {molecule_id}: {str(e)}")
            return []
            
    @classmethod
    def get_for_mixture(cls, mixture_id: str) -> Dict[str, List[Dict[str, Any]]]:
        """
        Get toxicity scores for components of a mixture.
        
        Args:
            mixture_id: UUID of the mixture
            
        Returns:
            Dictionary mapping molecule IDs to toxicity scores
        """
        try:
            # Get mixture components
            components_response = cls.get_supabase().table("mixture_component").select("molecule_id").eq("mixture_id", mixture_id).execute()
            
            if not components_response.data:
                return {}
                
            # Get toxicity scores for each component
            component_scores = {}
            for component in components_response.data:
                molecule_id = component["molecule_id"]
                scores = cls.get_for_molecule(molecule_id)
                if scores:
                    component_scores[molecule_id] = scores
                    
            return component_scores
            
        except Exception as e:
            logger.error(f"Error getting toxicity scores for mixture {mixture_id}: {str(e)}")
            return {}
            
    @classmethod
    def calculate_mixture_toxicity_score(cls, mixture_id: str) -> Dict[str, Any]:
        """
        Calculate aggregate toxicity score for a mixture based on component scores.
        
        Args:
            mixture_id: UUID of the mixture
            
        Returns:
            Dictionary with aggregate toxicity score
        """
        try:
            # Get mixture components with concentrations
            components_response = cls.get_supabase().table("mixture_component").select("*").eq("mixture_id", mixture_id).execute()
            
            if not components_response.data:
                return {"error": "No components found for this mixture"}
                
            # Get toxicity scores for each component
            component_scores = {}
            total_concentration = 0
            
            for component in components_response.data:
                molecule_id = component["molecule_id"]
                concentration = component["concentration"]
                total_concentration += concentration
                
                # Get overall toxicity score for the molecule
                molecule_response = cls.get_supabase().table("molecule").select("toxicity_score").eq("id", molecule_id).execute()
                
                if molecule_response.data and molecule_response.data[0].get("toxicity_score") is not None:
                    component_scores[molecule_id] = {
                        "toxicity_score": molecule_response.data[0]["toxicity_score"],
                        "concentration": concentration
                    }
            
            if not component_scores:
                return {"error": "No toxicity scores found for mixture components"}
                
            # Calculate weighted average based on concentration
            weighted_sum = 0
            for molecule_id, data in component_scores.items():
                weighted_sum += data["toxicity_score"] * (data["concentration"] / total_concentration)
                
            return {
                "mixture_id": mixture_id,
                "aggregate_toxicity_score": weighted_sum,
                "component_count": len(component_scores),
                "calculation_method": "concentration_weighted_average"
            }
            
        except Exception as e:
            logger.error(f"Error calculating mixture toxicity score for {mixture_id}: {str(e)}")
            return {"error": str(e)}

class UserProfile(BaseModel):
    """
    Model for user profiles in the database.
    Fields: id (uuid), user_id (Supabase UID), email, name, created_at, updated_at
    """
    table_name = "user_profile"

    @classmethod
    def create_or_update(cls, user_id: str, email: str, name: str = None) -> Dict[str, Any]:
        """
        Create or update a user profile for the given user_id.
        """
        supabase = get_supabase_client()
        now = datetime.utcnow().isoformat()
        data = {
            "user_id": user_id,
            "email": email,
            "name": name,
            "updated_at": now,
        }
        # Upsert by user_id
        response = supabase.table(cls.table_name).upsert(data, on_conflict="user_id").execute()
        if response.error:
            raise Exception(f"Error upserting user profile: {response.error}")
        return response.data[0] if response.data else {}

    @classmethod
    def get_by_user_id(cls, user_id: str) -> Optional[Dict[str, Any]]:
        supabase = get_supabase_client()
        response = supabase.table(cls.table_name).select("*").eq("user_id", user_id).single().execute()
        if response.error:
            return None
        return response.data

    @classmethod
    def update_profile(cls, user_id: str, updates: Dict[str, Any]) -> Dict[str, Any]:
        supabase = get_supabase_client()
        updates["updated_at"] = datetime.utcnow().isoformat()
        response = supabase.table(cls.table_name).update(updates).eq("user_id", user_id).execute()
        if response.error:
            raise Exception(f"Error updating user profile: {response.error}")
        return response.data[0] if response.data else {}

# Protocol-related fields
protocol_fields = {
    'id': fields.String,
    'mixture_id': fields.String,
    'name': fields.String,
    'description': fields.String,
    'target_concentration': fields.Float,
    'sample_type': fields.String,
    'starting_temperature': fields.Float,
    'target_temperature': fields.Float,
    'step_count': fields.Integer,
    'steps': fields.Raw,
    'custom_sensitivity': fields.Raw,
    'created_at': FlexibleDateTime(dt_format='iso8601'),
    'updated_at': FlexibleDateTime(dt_format='iso8601'),
    'created_by': fields.String,
    'mixture_name': fields.String
}

# Protocol-related schemas
class ProtocolSchema(Schema):
    """
    Schema for creating or updating a protocol.
    
    This schema validates the data for a cryopreservation protocol, including
    its basic information, temperature parameters, and procedural steps.
    """
    name = ma_fields.String(
        required=True,
        description="Name of the protocol",
        example="Slow Cooling Protocol for Stem Cells",
        validate=validate.Length(min=1, max=100)
    )
    description = ma_fields.String(
        description="Description of the protocol",
        example="Optimized cooling protocol for human mesenchymal stem cells",
        validate=validate.Length(max=500)
    )
    target_concentration = ma_fields.Float(
        required=True,
        validate=validate.Range(min=0),
        description="Target concentration of the cryoprotectant mixture (in appropriate units)",
        example=10.0
    )
    sample_type = ma_fields.String(
        required=True,
        description="Type of biological sample being preserved",
        example="human_mscs",
        validate=validate.Length(min=1, max=50)
    )
    starting_temperature = ma_fields.Float(
        required=True,
        description="Initial temperature in degrees Celsius",
        example=22.0
    )
    target_temperature = ma_fields.Float(
        description="Final target temperature in degrees Celsius",
        example=-196.0
    )
    step_count = ma_fields.Integer(
        validate=validate.Range(min=1),
        description="Number of steps in the protocol",
        example=5
    )
    steps = ma_fields.Raw(
        required=True,
        description="List of protocol steps with temperature, duration, and instructions",
        example=[
            {"temp": 22.0, "duration": 10, "instruction": "Add cryoprotectant solution"},
            {"temp": 4.0, "duration": 30, "instruction": "Cool to 4°C and hold"},
            {"temp": -80.0, "duration": 120, "instruction": "Transfer to -80°C freezer"}
        ]
    )
    custom_sensitivity = ma_fields.Raw(
        description="Custom sensitivity parameters for specific sample types",
        example={"temperature_threshold": -5.0, "cooling_rate_max": 1.0}
    )


class Protocol(BaseModel):
    """Model for protocols in the database."""
    
    table_name = 'protocols'
    
    @classmethod
    def create_protocol(cls, mixture_id: str, name: str, description: str,
                       target_concentration: float, sample_type: str,
                       starting_temperature: float, target_temperature: float = None,
                       step_count: int = None, steps: List[Dict[str, Any]] = None,
                       custom_sensitivity: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Create a new protocol.
        
        Args:
            mixture_id: ID of the mixture
            name: Name of the protocol
            description: Description of the protocol
            target_concentration: Target concentration
            sample_type: Type of biological sample
            starting_temperature: Initial temperature (°C)
            target_temperature: Final temperature (°C)
            step_count: Number of steps
            steps: List of protocol steps
            custom_sensitivity: Sensitivity profile
            
        Returns:
            The created protocol record
        """
        # Get the current user ID
        user_id = get_user_id()
        
        # Check if mixture exists
        mixture = Mixture.get(mixture_id)
        if not mixture:
            raise ValueError(f"Mixture with ID {mixture_id} not found")
        
        # Prepare protocol data
        protocol_data = {
            "mixture_id": mixture_id,
            "name": name,
            "description": description,
            "target_concentration": target_concentration,
            "sample_type": sample_type,
            "starting_temperature": starting_temperature,
            "target_temperature": target_temperature,
            "step_count": step_count,
            "steps": steps,
            "custom_sensitivity": custom_sensitivity,
            "created_by": user_id
        }
        
        # Insert protocol
        response = cls.get_table().insert(protocol_data).execute()
        
        if response.error:
            raise Exception(f"Error creating protocol: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def get_protocol(cls, id_value: str) -> Dict[str, Any]:
        """
        Get a protocol with mixture information.
        
        Args:
            id_value: The ID of the protocol
            
        Returns:
            Protocol with mixture information
        """
        response = cls.get_supabase().table("protocol_with_mixture").select("*").eq("id", id_value).execute()
        
        if response.error:
            raise Exception(f"Error getting protocol: {response.error}")
        
        return response.data[0] if response.data else None
    
    @classmethod
    def get_protocols_for_mixture(cls, mixture_id: str) -> List[Dict[str, Any]]:
        """
        Get all protocols for a mixture.
        
        Args:
            mixture_id: ID of the mixture
            
        Returns:
            List of protocols
        """
        response = cls.get_table().select('*').eq('mixture_id', mixture_id).execute()
        
        if response.error:
            raise Exception(f"Error getting protocols for mixture: {response.error}")
        
        return response.data
    
    @classmethod
    def update_protocol(cls, id_value: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update a protocol.
        
        Args:
            id_value: The ID of the protocol to update
            data: Dictionary containing the protocol data to update
            
        Returns:
            The updated protocol record
        """
        response = cls.get_table().update(data).eq(cls.id_field, id_value).execute()
        
        if response.error:
            raise Exception(f"Error updating protocol: {response.error}")
        
        return cls.get_protocol(id_value)
    
    
    class ProtocolStep:
        """
        Model for protocol steps.
        
        This class represents a single step in a cryopreservation protocol,
        including temperature, concentration, duration, and instructions.
        """
        
        @staticmethod
        def validate_step(step_data: Dict[str, Any]) -> Dict[str, Any]:
            """
            Validate a protocol step.
            
            Args:
                step_data: Dictionary containing step data
                
            Returns:
                Validated step data
                
            Raises:
                ValueError: If validation fails
            """
            required_fields = ["step", "action", "from_concentration", "to_concentration", "temperature", "hold_time_min"]
            for field in required_fields:
                if field not in step_data:
                    raise ValueError(f"Missing required field: {field}")
                    
            # Validate numeric fields
            if step_data["from_concentration"] < 0:
                raise ValueError("From concentration must be non-negative")
                
            if step_data["to_concentration"] < 0:
                raise ValueError("To concentration must be non-negative")
                
            if step_data["hold_time_min"] <= 0:
                raise ValueError("Hold time must be positive")
                
            return step_data
        
        @staticmethod
        def calculate_step_parameters(
            from_concentration: float,
            to_concentration: float,
            temperature: float,
            sample_type: str = "cell_line",
            custom_sensitivity: Optional[Dict[str, Any]] = None
        ) -> Dict[str, Any]:
            """
            Calculate optimal parameters for a protocol step.
            
            Args:
                from_concentration: Starting concentration
                to_concentration: Target concentration
                temperature: Temperature for the step
                sample_type: Type of biological sample
                custom_sensitivity: Custom sensitivity parameters
                
            Returns:
                Dictionary with calculated step parameters
            """
            # Get sensitivity profile for sample type
            sensitivity_profiles = {
                "cell_line": {
                    "osmotic_tolerance": 0.3,
                    "time_per_step": 5
                },
                "primary_cells": {
                    "osmotic_tolerance": 0.2,
                    "time_per_step": 10
                },
                "tissue": {
                    "osmotic_tolerance": 0.15,
                    "time_per_step": 15
                },
                "organoid": {
                    "osmotic_tolerance": 0.18,
                    "time_per_step": 12
                },
                "embryo": {
                    "osmotic_tolerance": 0.1,
                    "time_per_step": 8
                }
            }
            
            profile = sensitivity_profiles.get(sample_type, sensitivity_profiles["cell_line"])
            
            # Override with custom sensitivity if provided
            if custom_sensitivity:
                profile.update(custom_sensitivity)
                
            # Calculate hold time based on concentration change and sensitivity
            concentration_change = abs(to_concentration - from_concentration)
            base_time = profile["time_per_step"]
            
            # Adjust time based on concentration change
            adjusted_time = base_time * (1 + concentration_change / 10)
            
            # Round to nearest minute
            hold_time = max(1, round(adjusted_time))
            
            return {
                "hold_time_min": hold_time,
                "notes": f"Equilibrate for {hold_time} minutes. Adjust time based on sample response."
            }
    
    
    class ProtocolComparison:
        """
        Model for protocol comparisons.
        
        This class provides methods for comparing multiple protocols and
        analyzing their differences and similarities.
        """
        
        @staticmethod
        def create_comparison(protocol_ids: List[str]) -> Dict[str, Any]:
            """
            Create a comparison of multiple protocols.
            
            Args:
                protocol_ids: List of protocol IDs to compare
                
            Returns:
                Dictionary with comparison results
                
            Raises:
                ValueError: If fewer than 2 protocol IDs are provided
            """
            if len(protocol_ids) < 2:
                raise ValueError("At least two protocol IDs are required for comparison")
                
            # Get protocols
            protocols = []
            for protocol_id in protocol_ids:
                protocol = Protocol.get_protocol(protocol_id)
                if not protocol:
                    raise ValueError(f"Protocol with ID {protocol_id} not found")
                protocols.append(protocol)
                
            # Create comparison
            comparison = {
                "id": str(uuid.uuid4()),
                "protocol_ids": protocol_ids,
                "protocols": protocols,
                "created_at": datetime.now().isoformat(),
                "parameter_comparison": {},
                "step_comparison": {},
                "summary": {},
                "recommendations": []
            }
            
            # Compare parameters
            param_keys = ["target_concentration", "sample_type", "starting_temperature", "target_temperature", "step_count"]
            for key in param_keys:
                values = [p.get(key) for p in protocols]
                comparison["parameter_comparison"][key] = {
                    "values": values,
                    "same": len(set(str(v) for v in values)) == 1
                }
                
            # Compare steps
            max_steps = max(len(p.get("steps", [])) for p in protocols)
            for i in range(max_steps):
                step_comparison = {"step_number": i + 1, "data": []}
                for p in protocols:
                    steps = p.get("steps", [])
                    step_data = steps[i] if i < len(steps) else None
                    step_comparison["data"].append(step_data)
                comparison["step_comparison"][f"step_{i+1}"] = step_comparison
                
            # Generate summary
            comparison["summary"] = {
                "protocol_count": len(protocols),
                "sample_types": list(set(p.get("sample_type") for p in protocols)),
                "concentration_range": [
                    min(p.get("target_concentration", 0) for p in protocols),
                    max(p.get("target_concentration", 0) for p in protocols)
                ],
                "step_count_range": [
                    min(len(p.get("steps", [])) for p in protocols),
                    max(len(p.get("steps", [])) for p in protocols)
                ]
            }
            
            # Generate recommendations
            if len(set(p.get("sample_type") for p in protocols)) > 1:
                comparison["recommendations"].append(
                    "Protocols are designed for different sample types. Consider using the protocol "
                    "specifically designed for your sample type."
                )
                
            return comparison
        
        @staticmethod
        def get_comparison(comparison_id: str) -> Dict[str, Any]:
            """
            Get a protocol comparison by ID.
            
            Args:
                comparison_id: ID of the comparison
                
            Returns:
                Comparison data or None if not found
            """
            # In a real implementation, this would retrieve the comparison from a database
            # For now, we'll return a stub
            return None
        
        @staticmethod
        def save_comparison(comparison_data: Dict[str, Any]) -> Dict[str, Any]:
            """
            Save a protocol comparison.
            
            Args:
                comparison_data: Comparison data to save
                
            Returns:
                Saved comparison data
            """
            # In a real implementation, this would save the comparison to a database
            # For now, we'll return the input data with an ID
            comparison_data["id"] = str(uuid.uuid4())
            comparison_data["created_at"] = datetime.now().isoformat()
            return comparison_data
    
    
    # Import team models
    from api.team_models import (
    team_fields, team_member_fields, shared_resource_fields,
    comment_fields, notification_fields, activity_log_fields,
    TeamSchema, TeamMemberSchema, SharedResourceSchema,
    CommentSchema, NotificationSchema,
    Team, TeamMember, SharedResource, Comment, Notification, ActivityLog
)

class LabVerification(BaseModel):
    """
    Model for lab verification data.

    This model handles the creation, retrieval, and status updates of lab verifications
    for experiments. It provides methods to record a new verification, fetch verification
    details for a given experiment, and update the verification status.

    Table: lab_verifications
    """

    table_name = 'lab_verifications'

    # Verification states
    PENDING = 'pending'
    VERIFIED = 'verified'
    REJECTED = 'rejected'

    @classmethod
    def record_verification(cls, experiment_id: str,
                            verification_status: str,
                            verifier: str,
                            equipment_used: str,
                            comments: str = None) -> dict:
        """
        Record verification for an experiment.

        Args:
            experiment_id (str): The ID of the experiment being verified.
            verification_status (str): The status of verification ('pending', 'verified', 'rejected').
            verifier (str): The user or staff ID of the verifier.
            equipment_used (str): Description or ID of equipment used.
            comments (str, optional): Additional comments.

        Returns:
            dict: The created verification record.

        Raises:
            ValueError: If required fields are missing or invalid.
            Exception: For database errors.
        """
        if verification_status not in [cls.PENDING, cls.VERIFIED, cls.REJECTED]:
            raise ValueError("Invalid verification status.")
        if not experiment_id or not verifier or not equipment_used:
            raise ValueError("experiment_id, verifier, and equipment_used are required.")

        data = {
            "experiment_id": experiment_id,
            "verification_status": verification_status,
            "verifier": verifier,
            "equipment_used": equipment_used,
            "comments": comments
        }
        try:
            response = cls.get_table().insert(data).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Recording LabVerification for experiment {experiment_id}"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error recording LabVerification: {error_message}")
            return response.data[0] if response.data else None
        except Exception as e:
            context = f"Recording LabVerification for experiment {experiment_id}"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error recording LabVerification: {str(e)}")

    @classmethod
    def get_verification(cls, experiment_id: str) -> dict:
        """
        Get verification for an experiment.

        Args:
            experiment_id (str): The ID of the experiment.

        Returns:
            dict: The verification record, or None if not found.

        Raises:
            Exception: For database errors.
        """
        if not experiment_id:
            raise ValueError("experiment_id is required.")
        try:
            response = cls.get_table().select("*").eq("experiment_id", experiment_id).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Getting LabVerification for experiment {experiment_id}"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error getting LabVerification: {error_message}")
            return response.data[0] if response.data else None
        except Exception as e:
            context = f"Getting LabVerification for experiment {experiment_id}"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error getting LabVerification: {str(e)}")

    @classmethod
    def update_verification_status(cls, verification_id: str,
                                  new_status: str,
                                  comments: str = None) -> dict:
        """
        Update verification status.

        Args:
            verification_id (str): The ID of the verification record.
            new_status (str): The new status ('pending', 'verified', 'rejected').
            comments (str, optional): Additional comments.

        Returns:
            dict: The updated verification record.

        Raises:
            ValueError: If required fields are missing or invalid.
            Exception: For database errors.
        """
        if new_status not in [cls.PENDING, cls.VERIFIED, cls.REJECTED]:
            raise ValueError("Invalid new_status value.")
        if not verification_id:
            raise ValueError("verification_id is required.")

        update_data = {
            "verification_status": new_status
        }
        if comments is not None:
            update_data["comments"] = comments

        try:
            response = cls.get_table().update(update_data).eq("id", verification_id).execute()
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                context = f"Updating LabVerification status for id {verification_id}"
                handle_error(error_message, context=context, log_level='error')
                raise Exception(f"Error updating LabVerification status: {error_message}")
            return response.data[0] if response.data else None
        except Exception as e:
            context = f"Updating LabVerification status for id {verification_id}"
            handle_error(e, context=context, log_level='error')
            raise Exception(f"Error updating LabVerification status: {str(e)}")
