"""
CryoProtect Analyzer API - Data Models

This module contains data models for the API, including request and response schemas.
"""

from flask_restful import fields
from marshmallow import Schema, fields as ma_fields, validate, ValidationError, post_load

# Flask-RESTful response fields
molecule_fields = {
    'id': fields.String,
    'cid': fields.Integer,
    'name': fields.String,
    'molecular_formula': fields.String,
    'smiles': fields.String,
    'pubchem_link': fields.String,
    'created_at': fields.DateTime(dt_format='iso8601'),
    'updated_at': fields.DateTime(dt_format='iso8601'),
    'properties': fields.Raw
}

mixture_fields = {
    'id': fields.String,
    'name': fields.String,
    'description': fields.String,
    'created_at': fields.DateTime(dt_format='iso8601'),
    'updated_at': fields.DateTime(dt_format='iso8601'),
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
    'created_at': fields.DateTime(dt_format='iso8601')
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
    'created_at': fields.DateTime(dt_format='iso8601')
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

# Marshmallow validation schemas
class MoleculeComponentSchema(Schema):
    """Schema for a molecule component in a mixture."""
    molecule_id = ma_fields.UUID(required=True)
    concentration = ma_fields.Float(required=True, validate=validate.Range(min=0))
    concentration_unit = ma_fields.String(required=True)

class MixtureSchema(Schema):
    """Schema for creating or updating a mixture."""
    name = ma_fields.String(required=True)
    description = ma_fields.String()
    components = ma_fields.List(ma_fields.Nested(MoleculeComponentSchema), required=True)

    @post_load
    def validate_components(self, data, **kwargs):
        """Validate that the mixture has at least one component."""
        if not data.get('components'):
            raise ValidationError('Mixture must have at least one component')
        return data

class PropertyValueSchema(Schema):
    """Schema for a property value."""
    property_name = ma_fields.String(required=True)
    value = ma_fields.Raw(required=True)
    data_type = ma_fields.String(validate=validate.OneOf(['numeric', 'text', 'boolean']))

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
    """Schema for creating a prediction."""
    property_name = ma_fields.String(required=True)
    value = ma_fields.Raw(required=True)
    confidence = ma_fields.Float(required=True, validate=validate.Range(min=0, max=1))
    calculation_method = ma_fields.String(required=True)

class ExperimentSchema(Schema):
    """Schema for recording an experiment."""
    property_name = ma_fields.String(required=True)
    value = ma_fields.Raw(required=True)
    experimental_conditions = ma_fields.String()
    date_performed = ma_fields.Date(required=True)

class ComparisonQuerySchema(Schema):
    """Schema for querying a comparison."""
    property_name = ma_fields.String(required=True)

# RDKit-related schemas
class MoleculeInputSchema(Schema):
    """Schema for molecule input data."""
    molecule_data = ma_fields.String(required=True)
    input_format = ma_fields.String(validate=validate.OneOf(['smiles', 'mol', 'sdf']), default='smiles')

class MoleculeVisualizationSchema(Schema):
    """Schema for molecule visualization options."""
    molecule_data = ma_fields.String(required=True)
    input_format = ma_fields.String(validate=validate.OneOf(['smiles', 'mol', 'sdf']), default='smiles')
    width = ma_fields.Integer(default=400)
    height = ma_fields.Integer(default=300)
    highlight_atoms = ma_fields.List(ma_fields.Integer())

class SubstructureSearchSchema(Schema):
    """Schema for substructure search."""
    query_mol_data = ma_fields.String(required=True)
    target_mol_data = ma_fields.String(required=True)
    query_format = ma_fields.String(validate=validate.OneOf(['smarts', 'smiles']), default='smarts')
    target_format = ma_fields.String(validate=validate.OneOf(['smiles', 'mol', 'sdf']), default='smiles')

class SimilaritySearchSchema(Schema):
    """Schema for similarity search."""
    mol1_data = ma_fields.String(required=True)
    mol2_data = ma_fields.String(required=True)
    mol1_format = ma_fields.String(validate=validate.OneOf(['smiles', 'mol', 'sdf']), default='smiles')
    mol2_format = ma_fields.String(validate=validate.OneOf(['smiles', 'mol', 'sdf']), default='smiles')
    fingerprint_type = ma_fields.String(validate=validate.OneOf(['morgan', 'maccs', 'topological']), default='morgan')