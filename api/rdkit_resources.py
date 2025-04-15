"""
CryoProtect Analyzer API - RDKit Resources

This module contains API resources for molecular operations using RDKit.
"""

from flask import request, g, current_app
from flask_restful import Resource, marshal_with, abort

from api.models import (
    molecular_property_fields, molecule_visualization_fields,
    substructure_search_fields, similarity_fields,
    MoleculeInputSchema, MoleculeVisualizationSchema,
    SubstructureSearchSchema, SimilaritySearchSchema
)
from api.utils import token_required
from api.rdkit_utils import (
    calculate_all_properties, generate_molecule_svg,
    perform_substructure_search, calculate_similarity
)
from marshmallow import ValidationError


class MoleculePropertyResource(Resource):
    """Resource for calculating molecular properties using RDKit."""
    
    @marshal_with(molecular_property_fields)
    def post(self):
        """Calculate properties for a molecule."""
        try:
            # Validate request data
            schema = MoleculeInputSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Calculate properties
            properties = calculate_all_properties(
                data["molecule_data"],
                data.get("input_format", "smiles")
            )
            
            if "error" in properties:
                abort(400, message=properties["error"])
            
            return properties, 200
        except Exception as e:
            current_app.logger.error(f"Error calculating molecular properties: {str(e)}")
            abort(500, message=f"Error calculating molecular properties: {str(e)}")


class MoleculeVisualizationResource(Resource):
    """Resource for generating molecular visualizations."""
    
    @marshal_with(molecule_visualization_fields)
    def post(self):
        """Generate a visualization for a molecule."""
        try:
            # Validate request data
            schema = MoleculeVisualizationSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Generate SVG
            svg = generate_molecule_svg(
                data["molecule_data"],
                data.get("input_format", "smiles"),
                data.get("width", 400),
                data.get("height", 300),
                data.get("highlight_atoms")
            )
            
            if not svg:
                abort(400, message="Failed to generate molecule visualization")
            
            return {
                "svg": svg,
                "width": data.get("width", 400),
                "height": data.get("height", 300)
            }, 200
        except Exception as e:
            current_app.logger.error(f"Error generating molecule visualization: {str(e)}")
            abort(500, message=f"Error generating molecule visualization: {str(e)}")


class SubstructureSearchResource(Resource):
    """Resource for performing substructure searches."""
    
    @marshal_with(substructure_search_fields)
    def post(self):
        """Perform a substructure search."""
        try:
            # Validate request data
            schema = SubstructureSearchSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Perform search
            result = perform_substructure_search(
                data["query_mol_data"],
                data["target_mol_data"],
                data.get("query_format", "smarts"),
                data.get("target_format", "smiles")
            )
            
            if "error" in result:
                abort(400, message=result["error"])
            
            return result, 200
        except Exception as e:
            current_app.logger.error(f"Error performing substructure search: {str(e)}")
            abort(500, message=f"Error performing substructure search: {str(e)}")


class SimilaritySearchResource(Resource):
    """Resource for calculating molecular similarity."""
    
    @marshal_with(similarity_fields)
    def post(self):
        """Calculate similarity between two molecules."""
        try:
            # Validate request data
            schema = SimilaritySearchSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Calculate similarity
            result = calculate_similarity(
                data["mol1_data"],
                data["mol2_data"],
                data.get("mol1_format", "smiles"),
                data.get("mol2_format", "smiles"),
                data.get("fingerprint_type", "morgan")
            )
            
            if "error" in result:
                abort(400, message=result["error"])
            
            return result, 200
        except Exception as e:
            current_app.logger.error(f"Error calculating molecular similarity: {str(e)}")
            abort(500, message=f"Error calculating molecular similarity: {str(e)}")


class MoleculePropertyCalculationResource(Resource):
    """Resource for calculating properties for a specific molecule in the database."""
    
    @token_required
    def post(self, molecule_id):
        """Calculate and store properties for a molecule in the database."""
        try:
            from api.utils import get_supabase_client, get_user_id
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Get the molecule from the database
            response = supabase.table("molecules").select("*").eq("id", molecule_id).execute()
            
            if not response.data:
                abort(404, message=f"Molecule with ID {molecule_id} not found")
            
            molecule = response.data[0]
            smiles = molecule.get("smiles")
            
            if not smiles:
                abort(400, message="Molecule does not have SMILES data")
            
            # Calculate properties
            properties = calculate_all_properties(smiles, "smiles")
            
            if "error" in properties:
                abort(400, message=properties["error"])
            
            # Get property types
            response = supabase.table("property_types").select("id, name, data_type").execute()
            property_types = response.data
            
            # Prepare property inserts
            property_inserts = []
            
            # Flatten the properties dictionary for storage
            flattened_properties = {}
            
            # Add hydrogen bonding properties
            if "hydrogen_bonding" in properties:
                hb = properties["hydrogen_bonding"]
                flattened_properties["H-Bond Donors"] = hb.get("donors")
                flattened_properties["H-Bond Acceptors"] = hb.get("acceptors")
                flattened_properties["Total H-Bonds"] = hb.get("total")
            
            # Add other properties
            flattened_properties["LogP"] = properties.get("logp")
            flattened_properties["TPSA"] = properties.get("tpsa")
            
            # Add molecular properties
            if "molecular_properties" in properties:
                mp = properties["molecular_properties"]
                for key, value in mp.items():
                    flattened_properties[key.replace("_", " ").title()] = value
            
            # Add permeability properties
            if "permeability" in properties:
                perm = properties["permeability"]
                for key, value in perm.items():
                    flattened_properties[key.replace("_", " ").title()] = value
            
            # Create property inserts
            for property_name, value in flattened_properties.items():
                if value is None:
                    continue
                    
                property_type = next((pt for pt in property_types if pt["name"] == property_name), None)
                if not property_type:
                    # Skip properties that don't have a corresponding property type
                    continue
                
                property_insert = {
                    "molecule_id": molecule_id,
                    "property_type_id": property_type["id"],
                    "created_by": user_id
                }
                
                # Set the appropriate value field based on data type
                if property_type["data_type"] == "numeric":
                    try:
                        property_insert["numeric_value"] = float(value) if value is not None else None
                    except (ValueError, TypeError):
                        continue
                elif property_type["data_type"] == "text":
                    property_insert["text_value"] = str(value) if value is not None else None
                elif property_type["data_type"] == "boolean":
                    property_insert["boolean_value"] = bool(value) if value is not None else None
                
                property_inserts.append(property_insert)
            
            # Insert properties if we have any
            if property_inserts:
                response = supabase.table("molecular_properties").insert(property_inserts).execute()
                
                if response.error:
                    abort(400, message=f"Error inserting properties: {response.error}")
                
                return {"message": f"Calculated and stored {len(property_inserts)} properties for molecule {molecule_id}"}, 200
            else:
                return {"message": "No properties to store"}, 200
            
        except Exception as e:
            current_app.logger.error(f"Error calculating properties for molecule: {str(e)}")
            abort(500, message=f"Error calculating properties for molecule: {str(e)}")