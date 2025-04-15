"""
CryoProtect Analyzer API - Resources

This module contains the API resource classes that handle HTTP requests.
"""

from flask import request, g, current_app
from flask_restful import Resource, marshal_with, reqparse, abort
from marshmallow import ValidationError

from api.models import (
    molecule_fields, mixture_fields, prediction_fields, experiment_fields, comparison_fields,
    MixtureSchema, PredictionSchema, ExperimentSchema, ComparisonQuerySchema
)
from api.utils import get_supabase_client, token_required, handle_supabase_error, get_user_id


class MoleculeListResource(Resource):
    """Resource for listing and creating molecules."""
    
    @marshal_with(molecule_fields)
    def get(self):
        """Get a list of molecules with their properties."""
        try:
            supabase = get_supabase_client()
            response = supabase.table("molecule_with_properties").select("*").execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            return response.data, 200
        except Exception as e:
            current_app.logger.error(f"Error fetching molecules: {str(e)}")
            abort(500, message=f"Error fetching molecules: {str(e)}")
    
    @token_required
    @marshal_with(molecule_fields)
    def post(self):
        """Import a molecule from PubChem."""
        parser = reqparse.RequestParser()
        parser.add_argument('cid', type=int, required=True, help='PubChem Compound ID is required')
        args = parser.parse_args()
        
        try:
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Call the database function to import the molecule
            response = supabase.rpc(
                "import_molecule_from_pubchem",
                {"p_cid": args['cid'], "p_user_id": user_id}
            ).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            molecule_id = response.data
            
            # Get the imported molecule with properties
            response = supabase.table("molecule_with_properties").select("*").eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Molecule with ID {molecule_id} not found")
            
            return response.data[0], 201
        except Exception as e:
            current_app.logger.error(f"Error importing molecule: {str(e)}")
            abort(500, message=f"Error importing molecule: {str(e)}")


class MoleculeResource(Resource):
    """Resource for retrieving a specific molecule."""
    
    @marshal_with(molecule_fields)
    def get(self, molecule_id):
        """Get a molecule with its properties."""
        try:
            supabase = get_supabase_client()
            response = supabase.table("molecule_with_properties").select("*").eq("id", molecule_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Molecule with ID {molecule_id} not found")
            
            return response.data[0], 200
        except Exception as e:
            current_app.logger.error(f"Error fetching molecule: {str(e)}")
            abort(500, message=f"Error fetching molecule: {str(e)}")


class MixtureListResource(Resource):
    """Resource for listing and creating mixtures."""
    
    @marshal_with(mixture_fields)
    def get(self):
        """Get a list of mixtures with their components."""
        try:
            supabase = get_supabase_client()
            response = supabase.table("mixture_with_components").select("*").execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            return response.data, 200
        except Exception as e:
            current_app.logger.error(f"Error fetching mixtures: {str(e)}")
            abort(500, message=f"Error fetching mixtures: {str(e)}")
    
    @token_required
    @marshal_with(mixture_fields)
    def post(self):
        """Create a new mixture."""
        try:
            # Validate request data
            schema = MixtureSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Create mixture
            response = supabase.table("mixtures").insert({
                "name": data["name"],
                "description": data.get("description", ""),
                "created_by": user_id
            }).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            mixture_id = response.data[0]["id"]
            
            # Add components
            component_inserts = [{
                "mixture_id": mixture_id,
                "molecule_id": comp["molecule_id"],
                "concentration": comp["concentration"],
                "concentration_unit": comp["concentration_unit"],
                "created_by": user_id
            } for comp in data["components"]]
            
            response = supabase.table("mixture_components").insert(component_inserts).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                # Rollback mixture creation
                supabase.table("mixtures").delete().eq("id", mixture_id).execute()
                abort(status_code, message=error_message)
            
            # Get the created mixture with components
            response = supabase.table("mixture_with_components").select("*").eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Mixture with ID {mixture_id} not found")
            
            return response.data[0], 201
        except Exception as e:
            current_app.logger.error(f"Error creating mixture: {str(e)}")
            abort(500, message=f"Error creating mixture: {str(e)}")


class MixtureResource(Resource):
    """Resource for retrieving a specific mixture."""
    
    @marshal_with(mixture_fields)
    def get(self, mixture_id):
        """Get a mixture with its components."""
        try:
            supabase = get_supabase_client()
            response = supabase.table("mixture_with_components").select("*").eq("id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Mixture with ID {mixture_id} not found")
            
            return response.data[0], 200
        except Exception as e:
            current_app.logger.error(f"Error fetching mixture: {str(e)}")
            abort(500, message=f"Error fetching mixture: {str(e)}")


class PredictionResource(Resource):
    """Resource for creating and retrieving predictions."""
    
    @marshal_with(prediction_fields)
    def get(self, mixture_id):
        """Get predictions for a mixture."""
        try:
            supabase = get_supabase_client()
            
            # Join with property_types and calculation_methods to get names
            response = supabase.from_("predictions").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value, 
                confidence, created_at, 
                property_types(name), calculation_methods(name)
            """).eq("mixture_id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Transform the response to match the prediction_fields structure
            predictions = []
            for pred in response.data:
                prediction = {
                    'id': pred['id'],
                    'mixture_id': pred['mixture_id'],
                    'property_type_id': pred['property_type_id'],
                    'property_name': pred['property_types']['name'],
                    'numeric_value': pred['numeric_value'],
                    'text_value': pred['text_value'],
                    'boolean_value': pred['boolean_value'],
                    'confidence': pred['confidence'],
                    'calculation_method': pred['calculation_methods']['name'],
                    'created_at': pred['created_at']
                }
                predictions.append(prediction)
            
            return predictions, 200
        except Exception as e:
            current_app.logger.error(f"Error fetching predictions: {str(e)}")
            abort(500, message=f"Error fetching predictions: {str(e)}")
    
    @token_required
    @marshal_with(prediction_fields)
    def post(self, mixture_id):
        """Add a prediction for a mixture."""
        try:
            # Validate request data
            schema = PredictionSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Get property type ID
            response = supabase.table("property_types").select("id, data_type").eq("name", data["property_name"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Property type '{data['property_name']}' not found")
            
            property_type = response.data[0]
            
            # Get calculation method ID
            response = supabase.table("calculation_methods").select("id").eq("name", data["calculation_method"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Calculation method '{data['calculation_method']}' not found")
            
            calculation_method = response.data[0]
            
            # Prepare prediction insert
            prediction = {
                "mixture_id": mixture_id,
                "property_type_id": property_type["id"],
                "calculation_method_id": calculation_method["id"],
                "confidence": data["confidence"],
                "created_by": user_id
            }
            
            # Set the appropriate value field based on data type
            if property_type["data_type"] == "numeric":
                prediction["numeric_value"] = float(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "text":
                prediction["text_value"] = str(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "boolean":
                prediction["boolean_value"] = bool(data["value"]) if data["value"] is not None else None
            else:
                abort(400, message=f"Unknown data type '{property_type['data_type']}'")
            
            # Insert prediction
            response = supabase.table("predictions").upsert(
                prediction,
                on_conflict="mixture_id,property_type_id,calculation_method_id"
            ).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Get the created prediction with property and method names
            prediction_id = response.data[0]["id"]
            response = supabase.from_("predictions").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value, 
                confidence, created_at, 
                property_types(name), calculation_methods(name)
            """).eq("id", prediction_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Prediction with ID {prediction_id} not found")
            
            # Transform the response to match the prediction_fields structure
            pred = response.data[0]
            prediction = {
                'id': pred['id'],
                'mixture_id': pred['mixture_id'],
                'property_type_id': pred['property_type_id'],
                'property_name': pred['property_types']['name'],
                'numeric_value': pred['numeric_value'],
                'text_value': pred['text_value'],
                'boolean_value': pred['boolean_value'],
                'confidence': pred['confidence'],
                'calculation_method': pred['calculation_methods']['name'],
                'created_at': pred['created_at']
            }
            
            return prediction, 201
        except Exception as e:
            current_app.logger.error(f"Error adding prediction: {str(e)}")
            abort(500, message=f"Error adding prediction: {str(e)}")


class ExperimentResource(Resource):
    """Resource for creating and retrieving experiments."""
    
    @marshal_with(experiment_fields)
    def get(self, mixture_id):
        """Get experiments for a mixture."""
        try:
            supabase = get_supabase_client()
            
            # Join with property_types to get names
            response = supabase.from_("experiments").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value, 
                experimental_conditions, date_performed, created_at, 
                property_types(name)
            """).eq("mixture_id", mixture_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Transform the response to match the experiment_fields structure
            experiments = []
            for exp in response.data:
                experiment = {
                    'id': exp['id'],
                    'mixture_id': exp['mixture_id'],
                    'property_type_id': exp['property_type_id'],
                    'property_name': exp['property_types']['name'],
                    'numeric_value': exp['numeric_value'],
                    'text_value': exp['text_value'],
                    'boolean_value': exp['boolean_value'],
                    'experimental_conditions': exp['experimental_conditions'],
                    'date_performed': exp['date_performed'],
                    'created_at': exp['created_at']
                }
                experiments.append(experiment)
            
            return experiments, 200
        except Exception as e:
            current_app.logger.error(f"Error fetching experiments: {str(e)}")
            abort(500, message=f"Error fetching experiments: {str(e)}")
    
    @token_required
    @marshal_with(experiment_fields)
    def post(self, mixture_id):
        """Record an experiment for a mixture."""
        try:
            # Validate request data
            schema = ExperimentSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Get property type ID
            response = supabase.table("property_types").select("id, data_type").eq("name", data["property_name"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Property type '{data['property_name']}' not found")
            
            property_type = response.data[0]
            
            # Prepare experiment insert
            experiment = {
                "mixture_id": mixture_id,
                "property_type_id": property_type["id"],
                "experimental_conditions": data.get("experimental_conditions", ""),
                "date_performed": data["date_performed"].isoformat(),
                "created_by": user_id
            }
            
            # Set the appropriate value field based on data type
            if property_type["data_type"] == "numeric":
                experiment["numeric_value"] = float(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "text":
                experiment["text_value"] = str(data["value"]) if data["value"] is not None else None
            elif property_type["data_type"] == "boolean":
                experiment["boolean_value"] = bool(data["value"]) if data["value"] is not None else None
            else:
                abort(400, message=f"Unknown data type '{property_type['data_type']}'")
            
            # Insert experiment
            response = supabase.table("experiments").insert(experiment).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Get the created experiment with property name
            experiment_id = response.data[0]["id"]
            response = supabase.from_("experiments").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value, 
                experimental_conditions, date_performed, created_at, 
                property_types(name)
            """).eq("id", experiment_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Experiment with ID {experiment_id} not found")
            
            # Transform the response to match the experiment_fields structure
            exp = response.data[0]
            experiment = {
                'id': exp['id'],
                'mixture_id': exp['mixture_id'],
                'property_type_id': exp['property_type_id'],
                'property_name': exp['property_types']['name'],
                'numeric_value': exp['numeric_value'],
                'text_value': exp['text_value'],
                'boolean_value': exp['boolean_value'],
                'experimental_conditions': exp['experimental_conditions'],
                'date_performed': exp['date_performed'],
                'created_at': exp['created_at']
            }
            
            return experiment, 201
        except Exception as e:
            current_app.logger.error(f"Error recording experiment: {str(e)}")
            abort(500, message=f"Error recording experiment: {str(e)}")


class ComparisonResource(Resource):
    """Resource for comparing predictions with experiments."""
    
    @marshal_with(comparison_fields)
    def get(self, mixture_id):
        """Compare prediction with experiment for a mixture."""
        try:
            # Validate query parameters
            schema = ComparisonQuerySchema()
            try:
                args = schema.load(request.args)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            supabase = get_supabase_client()
            
            # Get property type ID
            response = supabase.table("property_types").select("id").eq("name", args["property_name"]).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message=f"Property type '{args['property_name']}' not found")
            
            property_type = response.data[0]
            
            # Call the database function to compare
            response = supabase.rpc(
                "compare_prediction_with_experiment",
                {
                    "p_mixture_id": mixture_id,
                    "p_property_type_id": property_type["id"]
                }
            ).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            return response.data, 200
        except Exception as e:
            current_app.logger.error(f"Error comparing prediction with experiment: {str(e)}")
            abort(500, message=f"Error comparing prediction with experiment: {str(e)}")