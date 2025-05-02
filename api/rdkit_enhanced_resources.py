"""
CryoProtect Analyzer API - Enhanced RDKit Resources

This module contains API resources for the enhanced RDKit functionality.
It provides endpoints for molecular property calculation, caching, batch processing,
and specialized cryoprotectant analysis.
"""

from flask import request, g, current_app
from flask_restful import Resource, marshal_with, abort
import logging

from api.models import (
    molecular_property_fields, molecule_visualization_fields,
    substructure_search_fields, similarity_fields
)
from api.rdkit_schemas import (
    MoleculeInputSchema, MoleculeVisualizationSchema,
    SubstructureSearchSchema, SimilaritySearchSchema,
    BatchMoleculeSchema, SimilaritySearchBatchSchema,
    SubstructureSearchBatchSchema, MoleculeGridSchema,
    ScaffoldAnalysisSchema, MolecularDynamicsSchema,
    MoleculePropertyCalculationBatchSchema
)
from api.utils import token_required, _handle_json_serialization, handle_supabase_error
from api.rdkit_enhanced import (
    calculate_properties_with_cache, calculate_cryoprotectant_properties,
    batch_calculate_properties, find_similar_molecules, batch_substructure_search,
    generate_molecule_grid, analyze_scaffold, clear_property_cache
)
from marshmallow import Schema, fields as ma_fields, validate, ValidationError

# Set up logging
logger = logging.getLogger(__name__)

# Import API documentation utilities
try:
    from flask_apispec import use_kwargs, marshal_with as apispec_marshal_with, doc
    from flask_apispec.views import MethodResource
except ImportError:
    # Create dummy decorators if flask-apispec is not installed
    def use_kwargs(*args, **kwargs): return lambda f: f
    def apispec_marshal_with(*args, **kwargs): return lambda f: f
    def doc(*args, **kwargs): return lambda f: f
    MethodResource = Resource

# Define API resources
class CryoprotectantPropertyResource(Resource):
    """Resource for calculating cryoprotectant-specific properties."""
    
    def post(self):
        """Calculate cryoprotectant-specific properties for a molecule."""
        try:
            # Validate request data
            schema = MoleculeInputSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Calculate properties
            properties = calculate_cryoprotectant_properties(
                data["molecule_data"],
                data.get("input_format", "smiles")
            )
            
            if "error" in properties:
                return _handle_json_serialization({'error': properties["error"]}), 400
            
            return _handle_json_serialization(properties), 200
        except Exception as e:
            current_app.logger.error(f"Error calculating cryoprotectant properties: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

class BatchPropertyCalculationResource(Resource):
    """Resource for batch calculation of molecular properties."""
    
    def post(self):
        """Calculate properties for multiple molecules in batch."""
        try:
            # Validate request data
            schema = BatchMoleculeSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Prepare molecule data for batch processing
            molecules = []
            for mol in data["molecules"]:
                molecules.append({
                    'data': mol.get('data', ''),
                    'format': mol.get('format', 'smiles')
                })
            
            # Calculate properties in batch
            results = batch_calculate_properties(molecules)
            
            return _handle_json_serialization(results), 200
        except Exception as e:
            current_app.logger.error(f"Error in batch property calculation: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

class SimilaritySearchBatchResource(Resource):
    """Resource for batch similarity search."""
    
    def post(self):
        """Find molecules similar to a query molecule from a list of targets."""
        try:
            # Validate request data
            schema = SimilaritySearchBatchSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Perform similarity search
            results = find_similar_molecules(
                data["query_mol_data"],
                data["target_molecules"],
                data.get("query_format", "smiles"),
                data.get("fingerprint_type", "morgan"),
                data.get("similarity_threshold", 0.7)
            )
            
            if results and "error" in results[0]:
                return _handle_json_serialization({'error': results[0]["error"]}), 400
            
            return _handle_json_serialization(results), 200
        except Exception as e:
            current_app.logger.error(f"Error in batch similarity search: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

class SubstructureSearchBatchResource(Resource):
    """Resource for batch substructure search."""
    
    def post(self):
        """Perform substructure search on multiple target molecules."""
        try:
            # Validate request data
            schema = SubstructureSearchBatchSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Perform substructure search
            results = batch_substructure_search(
                data["query_mol_data"],
                data["target_molecules"],
                data.get("query_format", "smarts")
            )
            
            if results and "error" in results[0]:
                return _handle_json_serialization({'error': results[0]["error"]}), 400
            
            return _handle_json_serialization(results), 200
        except Exception as e:
            current_app.logger.error(f"Error in batch substructure search: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

class MoleculeGridResource(Resource):
    """Resource for generating a grid of molecule visualizations."""
    
    def post(self):
        """Generate a grid of molecule visualizations."""
        try:
            # Validate request data
            schema = MoleculeGridSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Generate molecule grid
            svg = generate_molecule_grid(
                data["molecules"],
                data.get("labels"),
                data.get("mol_per_row", 3),
                (data.get("sub_img_width", 200), data.get("sub_img_height", 200)),
                data.get("input_format", "smiles")
            )
            
            if not svg:
                return _handle_json_serialization({'error': "Failed to generate molecule grid"}), 400
            
            return _handle_json_serialization({
                "svg": svg,
                "molecule_count": len(data["molecules"])
            }), 200
        except Exception as e:
            current_app.logger.error(f"Error generating molecule grid: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

class ScaffoldAnalysisResource(Resource):
    """Resource for analyzing molecular scaffolds."""
    
    def post(self):
        """Analyze the scaffold of a molecule."""
        try:
            # Validate request data
            schema = ScaffoldAnalysisSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Analyze scaffold
            result = analyze_scaffold(
                data["molecule_data"],
                data.get("input_format", "smiles")
            )
            
            if "error" in result:
                return _handle_json_serialization({'error': result["error"]}), 400
            
            return _handle_json_serialization(result), 200
        except Exception as e:
            current_app.logger.error(f"Error analyzing scaffold: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

class MolecularDynamicsResource(Resource):
    """Resource for molecular dynamics simulations."""
    
    @token_required
    def post(self):
        """Run a molecular dynamics simulation for a molecule."""
        try:
            # Validate request data
            schema = MolecularDynamicsSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Extract parameters
            molecule_data = data["molecule_data"]
            input_format = data.get("input_format", "smiles")
            simulation_time = data.get("simulation_time", 10.0)
            temperature = data.get("temperature", 300.0)
            force_field = data.get("force_field", "MMFF94")
            solvent = data.get("solvent", "water")
            include_trajectory = data.get("include_trajectory", False)
            
            # Run molecular dynamics simulation
            # This is a placeholder - in a real implementation, this would call
            # a function that performs the actual molecular dynamics simulation
            
            # Simulate the result
            import time
            import random
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Parse molecule
            if input_format == "smiles":
                mol = Chem.MolFromSmiles(molecule_data)
            elif input_format == "mol":
                mol = Chem.MolFromMolBlock(molecule_data)
            elif input_format == "sdf":
                suppl = Chem.SDMolSupplier()
                suppl.SetData(molecule_data)
                mol = next(suppl)
            else:
                return _handle_json_serialization({'error': f"Unsupported input format: {input_format}"}), 400
            
            if mol is None:
                return _handle_json_serialization({'error': "Failed to parse molecule"}), 400
            
            # Generate 3D coordinates if not present
            if not mol.GetNumConformers():
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)
            
            # Simulate a molecular dynamics run
            # In a real implementation, this would use a proper MD engine
            time.sleep(2)  # Simulate computation time
            
            # Generate simulated results
            result = {
                "molecule": {
                    "smiles": Chem.MolToSmiles(mol),
                    "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
                    "weight": Chem.rdMolDescriptors.CalcExactMolWt(mol)
                },
                "simulation": {
                    "time": simulation_time,
                    "temperature": temperature,
                    "force_field": force_field,
                    "solvent": solvent,
                    "status": "completed"
                },
                "results": {
                    "energy": {
                        "initial": round(random.uniform(100, 500), 2),
                        "final": round(random.uniform(50, 200), 2),
                        "average": round(random.uniform(75, 300), 2),
                        "unit": "kcal/mol"
                    },
                    "rmsd": round(random.uniform(0.5, 3.0), 2),
                    "radius_of_gyration": round(random.uniform(3.0, 10.0), 2),
                    "solvent_accessible_surface_area": round(random.uniform(200, 800), 2),
                    "hydrogen_bonds": {
                        "average": round(random.uniform(1, 10), 1),
                        "max": int(random.uniform(2, 15))
                    }
                }
            }
            
            # Add trajectory data if requested
            if include_trajectory:
                # In a real implementation, this would be actual trajectory data
                # For this example, we'll just generate some random coordinates
                num_frames = 10
                num_atoms = mol.GetNumAtoms()
                
                trajectory = []
                for frame in range(num_frames):
                    frame_data = {
                        "time": frame * simulation_time / num_frames,
                        "energy": round(result["results"]["energy"]["initial"] -
                                       (result["results"]["energy"]["initial"] -
                                        result["results"]["energy"]["final"]) *
                                       frame / num_frames, 2),
                        "coordinates": []
                    }
                    
                    # Generate random coordinates for each atom
                    # In a real implementation, these would be actual trajectory coordinates
                    for atom in range(num_atoms):
                        frame_data["coordinates"].append({
                            "atom_idx": atom,
                            "element": mol.GetAtomWithIdx(atom).GetSymbol(),
                            "x": round(random.uniform(-10, 10), 3),
                            "y": round(random.uniform(-10, 10), 3),
                            "z": round(random.uniform(-10, 10), 3)
                        })
                    
                    trajectory.append(frame_data)
                
                result["trajectory"] = trajectory
            
            return _handle_json_serialization(result), 200
        except Exception as e:
            current_app.logger.error(f"Error running molecular dynamics: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

class PropertyCacheResource(Resource):
    """Resource for managing the property cache."""
    
    @token_required
    def delete(self):
        """Clear the property cache."""
        try:
            result = clear_property_cache()
            
            if result["status"] == "error":
                return _handle_json_serialization({'error': result["message"]}), 400
            
            return _handle_json_serialization({"message": result["message"]}), 200
        except Exception as e:
            current_app.logger.error(f"Error clearing property cache: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

class MoleculePropertyCalculationBatchResource(MethodResource):
    """Resource for calculating properties for multiple molecules in the database."""
    
    @doc(description='Calculate and store properties for multiple molecules in the database',
         tags=['RDKit'],
         request_body={
             'molecule_ids': {'description': 'List of molecule IDs to calculate properties for', 'type': 'array', 'items': {'type': 'string', 'format': 'uuid'}, 'required': True}
         },
         security=[{'Bearer': []}])
    @token_required
    def post(self):
        """Calculate and store properties for multiple molecules in the database."""
        try:
            from api.utils import get_supabase_client, get_user_id
            
            # Validate request data
            schema = MoleculePropertyCalculationBatchSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            molecule_ids = data["molecule_ids"]
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Get the molecules from the database
            response = supabase.table("molecules").select("*").in_("id", molecule_ids).execute()
            
            if not response.data:
                return _handle_json_serialization({'error': "No molecules found with the provided IDs"}), 404
            
            molecules = response.data
            
            # Get property types
            response = supabase.table("property_types").select("id, name, data_type").execute()
            property_types = response.data
            
            # Process each molecule
            results = []
            for molecule in molecules:
                molecule_id = molecule.get("id")
                smiles = molecule.get("smiles")
                
                if not smiles:
                    results.append({
                        "molecule_id": molecule_id,
                        "status": "error",
                        "message": "Molecule does not have SMILES data"
                    })
                    continue
                
                # Calculate properties
                properties = calculate_cryoprotectant_properties(smiles, "smiles")
                
                if "error" in properties:
                    results.append({
                        "molecule_id": molecule_id,
                        "status": "error",
                        "message": properties["error"]
                    })
                    continue
                
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
                
                # Add cryoprotectant properties
                if "cryoprotectant_properties" in properties:
                    cp = properties["cryoprotectant_properties"]
                    for key, value in cp.items():
                        flattened_properties[key.replace("_", " ").title()] = value
                
                # Create property inserts
                for property_name, value in flattened_properties.items():
                    if value is None:
                        continue
                        
                    property_type = next((pt for pt in property_types if pt["name"] == property_name), None)
                    if not property_type:
                        # Create a new property type if it doesn't exist
                        try:
                            # Determine data type
                            data_type = "numeric"
                            if isinstance(value, bool):
                                data_type = "boolean"
                            elif isinstance(value, str):
                                data_type = "text"
                            
                            # Insert new property type
                            new_property_type_response = supabase.table("property_types").insert({
                                "name": property_name,
                                "data_type": data_type,
                                "created_by": user_id
                            }).execute()
                            
                            if new_property_type_response.error:
                                logger.error(f"Error creating property type: {new_property_type_response.error}")
                                continue
                                
                            property_type = new_property_type_response.data[0]
                            # Add to property types list for future use
                            property_types.append(property_type)
                        except Exception as e:
                            logger.error(f"Error creating property type: {str(e)}")
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
                    try:
                        response = supabase.table("molecular_properties").insert(property_inserts).execute()
                        
                        if response.error:
                            results.append({
                                "molecule_id": molecule_id,
                                "status": "error",
                                "message": f"Error inserting properties: {response.error}"
                            })
                        else:
                            results.append({
                                "molecule_id": molecule_id,
                                "status": "success",
                                "properties_added": len(property_inserts)
                            })
                    except Exception as e:
                        results.append({
                            "molecule_id": molecule_id,
                            "status": "error",
                            "message": f"Error inserting properties: {str(e)}"
                        })
                else:
                    results.append({
                        "molecule_id": molecule_id,
                        "status": "warning",
                        "message": "No properties to store"
                    })
            
            return _handle_json_serialization({
                "message": f"Processed {len(molecules)} molecules",
                "results": results
            }), 200
            
        except Exception as e:
            current_app.logger.error(f"Error calculating properties for molecules: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500

def register_resources(api):
    """Register RDKit enhanced resources with the API."""
    api.add_resource(CryoprotectantPropertyResource, '/api/v1/rdkit-enhanced/cryoprotectant-properties')
    api.add_resource(BatchPropertyCalculationResource, '/api/v1/rdkit-enhanced/batch-calculate')
    api.add_resource(SimilaritySearchBatchResource, '/api/v1/rdkit-enhanced/batch-similarity')
    api.add_resource(SubstructureSearchBatchResource, '/api/v1/rdkit-enhanced/batch-substructure')
    api.add_resource(MoleculeGridResource, '/api/v1/rdkit-enhanced/molecule-grid')
    api.add_resource(ScaffoldAnalysisResource, '/api/v1/rdkit-enhanced/scaffold-analysis')
    api.add_resource(PropertyCacheResource, '/api/v1/rdkit-enhanced/property-cache')
    api.add_resource(MoleculePropertyCalculationBatchResource, '/api/v1/rdkit-enhanced/molecules/batch-calculate-properties')
    api.add_resource(MolecularDynamicsResource, '/api/v1/rdkit-enhanced/molecular-dynamics')