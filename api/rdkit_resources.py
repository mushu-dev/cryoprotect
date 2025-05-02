"""
CryoProtect Analyzer API - RDKit Resources

This module contains API resources for molecular operations using RDKit, including:
- Calculating molecular properties
- Generating molecular visualizations
- Performing substructure searches
- Calculating molecular similarity
- Storing calculated properties in the database

All endpoints follow RESTful design principles and return standardized
response formats. Authentication is required for operations that modify
the database.
"""

from flask import request, g, current_app
from flask_restful import Resource, marshal_with, abort
from typing import Dict, List, Any, Tuple, Optional, Union

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

from api.models import (
    molecular_property_fields, molecule_visualization_fields,
    substructure_search_fields, similarity_fields
)
from api.rdkit_schemas import (
    MoleculeInputSchema, MoleculeVisualizationSchema,
    SubstructureSearchSchema, SimilaritySearchSchema
)
from api.utils import token_required, _handle_json_serialization, handle_supabase_error
from api.rdkit_utils import (
    calculate_all_properties, generate_molecule_svg,
    perform_substructure_search, calculate_similarity
)
from marshmallow import ValidationError


class MoleculePropertyResource(MethodResource):
    """Resource for calculating molecular properties using RDKit.
    
    This endpoint calculates various physicochemical and structural properties
    of molecules using RDKit. Properties include logP, TPSA, hydrogen bonding,
    molecular weight, and other descriptors relevant for cryoprotection analysis.
    """
    
    @doc(description='Calculate comprehensive physicochemical and structural properties for a molecule using RDKit',
         tags=['RDKit'],
         request_body={
             'molecule_data': {'description': 'Molecule data (SMILES, InChI, etc.)', 'type': 'string', 'required': True, 'example': 'CC(=O)OC1=CC=CC=C1C(=O)O'},
             'input_format': {'description': 'Format of the input data', 'type': 'string', 'default': 'smiles', 'enum': ['smiles', 'inchi', 'mol'], 'example': 'smiles'}
         },
         responses={
             '200': {
                 'description': 'Calculated molecular properties',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'hydrogen_bonding': {
                                     'type': 'object',
                                     'properties': {
                                         'donors': {'type': 'integer', 'description': 'Number of hydrogen bond donors', 'example': 1},
                                         'acceptors': {'type': 'integer', 'description': 'Number of hydrogen bond acceptors', 'example': 4},
                                         'total': {'type': 'integer', 'description': 'Total hydrogen bond count', 'example': 5}
                                     }
                                 },
                                 'logp': {'type': 'number', 'description': 'Octanol-water partition coefficient (lipophilicity)', 'example': 1.43},
                                 'tpsa': {'type': 'number', 'description': 'Topological polar surface area in Å²', 'example': 63.6},
                                 'molecular_properties': {
                                     'type': 'object',
                                     'properties': {
                                         'molecular_weight': {'type': 'number', 'description': 'Molecular weight in g/mol', 'example': 180.16},
                                         'exact_mass': {'type': 'number', 'description': 'Exact mass in g/mol', 'example': 180.04},
                                         'heavy_atom_count': {'type': 'integer', 'description': 'Number of non-hydrogen atoms', 'example': 13},
                                         'atom_count': {'type': 'integer', 'description': 'Total number of atoms', 'example': 21},
                                         'rotatable_bond_count': {'type': 'integer', 'description': 'Number of rotatable bonds', 'example': 3},
                                         'ring_count': {'type': 'integer', 'description': 'Number of rings', 'example': 1},
                                         'aromatic_ring_count': {'type': 'integer', 'description': 'Number of aromatic rings', 'example': 1},
                                         'fraction_csp3': {'type': 'number', 'description': 'Fraction of sp3 hybridized carbons', 'example': 0.22}
                                     }
                                 },
                                 'functional_groups': {
                                     'type': 'object',
                                     'description': 'Counts of functional groups present in the molecule',
                                     'example': {'ester': 1, 'carboxylic_acid': 1, 'aromatic': 1}
                                 },
                                 'permeability': {
                                     'type': 'object',
                                     'properties': {
                                         'rule_of_5_violations': {'type': 'integer', 'description': 'Number of Lipinski Rule of 5 violations', 'example': 0},
                                         'veber_violations': {'type': 'integer', 'description': 'Number of Veber rule violations', 'example': 0},
                                         'bbb_permeant': {'type': 'boolean', 'description': 'Blood-brain barrier permeability prediction', 'example': True},
                                         'intestinal_absorption': {'type': 'boolean', 'description': 'Intestinal absorption prediction', 'example': True},
                                         'estimated_log_papp': {'type': 'number', 'description': 'Estimated log of apparent permeability', 'example': -4.1}
                                     }
                                 },
                                 'smiles': {'type': 'string', 'description': 'SMILES representation of the molecule', 'example': 'CC(=O)OC1=CC=CC=C1C(=O)O'},
                                 'inchi': {'type': 'string', 'description': 'InChI representation of the molecule', 'example': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'},
                                 'inchi_key': {'type': 'string', 'description': 'InChI key for the molecule', 'example': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'}
                             }
                         },
                         'examples': {
                             'Aspirin': {
                                 'value': {
                                     'hydrogen_bonding': {'donors': 1, 'acceptors': 4, 'total': 5},
                                     'logp': 1.43,
                                     'tpsa': 63.6,
                                     'molecular_properties': {
                                         'molecular_weight': 180.16,
                                         'exact_mass': 180.04,
                                         'heavy_atom_count': 13,
                                         'atom_count': 21,
                                         'rotatable_bond_count': 3,
                                         'ring_count': 1,
                                         'aromatic_ring_count': 1,
                                         'fraction_csp3': 0.22,
                                         'molecular_volume': 163.12
                                     },
                                     'functional_groups': {
                                         'ester': 1,
                                         'carboxylic_acid': 1,
                                         'aromatic': 1
                                     },
                                     'permeability': {
                                         'rule_of_5_violations': 0,
                                         'veber_violations': 0,
                                         'bbb_permeant': True,
                                         'intestinal_absorption': True,
                                         'estimated_log_papp': -4.1
                                     },
                                     'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                                     'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
                                     'inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
                                 }
                             },
                             'Glycerol': {
                                 'value': {
                                     'hydrogen_bonding': {'donors': 3, 'acceptors': 3, 'total': 6},
                                     'logp': -1.76,
                                     'tpsa': 60.69,
                                     'molecular_properties': {
                                         'molecular_weight': 92.09,
                                         'exact_mass': 92.05,
                                         'heavy_atom_count': 6,
                                         'atom_count': 14,
                                         'rotatable_bond_count': 2,
                                         'ring_count': 0,
                                         'aromatic_ring_count': 0,
                                         'fraction_csp3': 1.0,
                                         'molecular_volume': 88.12
                                     },
                                     'functional_groups': {
                                         'alcohol': 3,
                                         'hydroxyl': 3
                                     },
                                     'permeability': {
                                         'rule_of_5_violations': 0,
                                         'veber_violations': 0,
                                         'bbb_permeant': False,
                                         'intestinal_absorption': True,
                                         'estimated_log_papp': -5.9
                                     },
                                     'smiles': 'C(C(CO)O)O',
                                     'inchi': 'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
                                     'inchi_key': 'PEDCQBHIVMGVHV-UHFFFAOYSA-N'
                                 }
                             }
                         }
                     }
                 }
             },
             '400': {
                 'description': 'Invalid request data or molecule cannot be processed',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'error': {'type': 'string'}
                             }
                         },
                         'example': {
                             'error': 'Failed to parse molecule'
                         }
                     }
                 }
             },
             '500': {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'error': {'type': 'string'}
                             }
                         },
                         'example': {
                             'error': 'Internal server error'
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(molecular_property_fields, code=200, description='Calculated molecular properties')
    @marshal_with(molecular_property_fields)
    def post(self) -> Tuple[Dict[str, Any], int]:
        """Calculate comprehensive properties for a molecule.
        
        Calculates a comprehensive set of molecular properties using RDKit,
        including physical properties, structural features, and descriptors
        relevant for cryoprotection analysis. This endpoint is particularly
        useful for analyzing potential cryoprotectant molecules and predicting
        their behavior in solution.
        
        Returns:
            Tuple[Dict[str, Any], int]: Calculated properties and HTTP status code
            
        Raises:
            400: If request data is invalid or molecule cannot be processed
            500: If server error occurs
            
        Request Body:
            molecule_data (str): Molecule data in the specified format (required)
                Example: "CC(=O)OC1=CC=CC=C1C(=O)O" (Aspirin)
            input_format (str): Format of the input data (default: 'smiles')
                Supported formats: 'smiles', 'inchi', 'mol'
                
        Response:
            Properties include:
            - hydrogen_bonding: Dictionary containing:
                - donors: Number of hydrogen bond donors
                - acceptors: Number of hydrogen bond acceptors
                - total: Total hydrogen bond count
            - logp: Octanol-water partition coefficient (lipophilicity)
            - tpsa: Topological polar surface area in Å²
            - molecular_properties: Dictionary containing:
                - molecular_weight: Molecular weight in g/mol
                - exact_mass: Exact mass in g/mol
                - heavy_atom_count: Number of non-hydrogen atoms
                - atom_count: Total number of atoms
                - rotatable_bond_count: Number of rotatable bonds
                - ring_count: Number of rings
                - aromatic_ring_count: Number of aromatic rings
                - fraction_csp3: Fraction of sp3 hybridized carbons
                - molecular_volume: Molecular volume in Å³ (if 3D coordinates available)
            - functional_groups: Dictionary of functional groups present in the molecule
                with their counts (e.g., alcohol, amine, carboxylic_acid)
            - permeability: Dictionary containing:
                - rule_of_5_violations: Number of Lipinski Rule of 5 violations
                - veber_violations: Number of Veber rule violations
                - bbb_permeant: Blood-brain barrier permeability prediction (boolean)
                - intestinal_absorption: Intestinal absorption prediction (boolean)
                - estimated_log_papp: Estimated log of apparent permeability
            - smiles: SMILES representation of the molecule
            - inchi: InChI representation of the molecule
            - inchi_key: InChI key for the molecule
            
        Developer Notes:
            - This endpoint performs computationally intensive calculations and may take
              longer for complex molecules
            - For batch processing of multiple molecules, consider using a dedicated batch endpoint
            - The 3D molecular properties (like molecular_volume) require generating 3D coordinates,
              which may not be accurate for all molecules
            - All property calculations are performed using RDKit's implementation of standard
              cheminformatics algorithms
        """
        try:
            # Validate request data
            schema = MoleculeInputSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Calculate properties
            properties = calculate_all_properties(
                data["molecule_data"],
                data.get("input_format", "smiles")
            )
            
            if "error" in properties:
                return _handle_json_serialization({'error': properties["error"]}), 400
            
            return _handle_json_serialization(properties), 200
        except Exception as e:
            current_app.logger.error(f"Error calculating molecular properties: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500


class MoleculeVisualizationResource(MethodResource):
    """Resource for generating molecular visualizations.
    
    This endpoint generates SVG visualizations of molecules using RDKit.
    It supports customization of size and highlighting specific atoms
    to create informative and visually appealing molecular representations.
    """
    
    @doc(description='Generate an SVG visualization for a molecule with customizable dimensions and optional atom highlighting',
         tags=['RDKit'],
         request_body={
             'molecule_data': {'description': 'Molecule data (SMILES, InChI, etc.)', 'type': 'string', 'required': True, 'example': 'C1=CC=C(C=C1)C(=O)O'},
             'input_format': {'description': 'Format of the input data', 'type': 'string', 'default': 'smiles', 'enum': ['smiles', 'inchi', 'mol'], 'example': 'smiles'},
             'width': {'description': 'Width of the SVG in pixels', 'type': 'integer', 'default': 400, 'example': 400},
             'height': {'description': 'Height of the SVG in pixels', 'type': 'integer', 'default': 300, 'example': 300},
             'highlight_atoms': {'description': 'List of atom indices to highlight (0-based)', 'type': 'array', 'items': {'type': 'integer'}, 'example': [0, 1, 2]}
         },
         responses={
             '200': {
                 'description': 'Molecule visualization as SVG',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'svg': {'type': 'string', 'description': 'SVG string representation of the molecule'},
                                 'width': {'type': 'integer', 'description': 'Width of the SVG in pixels', 'example': 400},
                                 'height': {'type': 'integer', 'description': 'Height of the SVG in pixels', 'example': 300}
                             }
                         },
                         'examples': {
                             'Benzoic Acid': {
                                 'value': {
                                     'svg': '<?xml version="1.0" encoding="UTF-8"?>\n<svg baseProfile="full" height="300px" version="1.1" viewBox="0 0 400 300" width="400px" xml:space="preserve" xmlns="http://www.w3.org/2000/svg" xmlns:rdkit="http://www.rdkit.org/xml" xmlns:xlink="http://www.w3.org/1999/xlink">\n<!-- END OF HEADER -->\n<rect height="300" style="opacity:1.0;fill:none;stroke:none" width="400" x="0" y="0">\n</rect>\n<path d="M 183.636,147.879 L 183.636,147.879" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,147.879 L 206.364,132.121" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,147.879 L 160.909,132.121" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,147.879 L 183.636,179.394" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 206.364,132.121 L 229.091,147.879" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 160.909,132.121 L 138.182,147.879" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 229.091,147.879 L 229.091,179.394" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 138.182,147.879 L 138.182,179.394" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 229.091,179.394 L 206.364,195.152" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 138.182,179.394 L 160.909,195.152" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 206.364,195.152 L 183.636,179.394" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 160.909,195.152 L 183.636,179.394" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,179.394 L 183.636,210.909" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,210.909 L 206.364,226.667" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,210.909 L 160.909,226.667" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 206.364,226.667 L 206.364,226.667" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 206.364,226.667 L 229.091,210.909" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 160.909,226.667 L 160.909,226.667" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 160.909,226.667 L 138.182,210.909" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 229.091,210.909 L 229.091,210.909" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 138.182,210.909 L 138.182,210.909" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n</svg>',
                                     'width': 400,
                                     'height': 300
                                 }
                             },
                             'Aspirin with Highlighted Atoms': {
                                 'value': {
                                     'svg': '<?xml version="1.0" encoding="UTF-8"?>\n<svg baseProfile="full" height="300px" version="1.1" viewBox="0 0 400 300" width="400px" xml:space="preserve" xmlns="http://www.w3.org/2000/svg" xmlns:rdkit="http://www.rdkit.org/xml" xmlns:xlink="http://www.w3.org/1999/xlink">\n<!-- END OF HEADER -->\n<rect height="300" style="opacity:1.0;fill:none;stroke:none" width="400" x="0" y="0">\n</rect>\n<ellipse cx="183.636" cy="147.879" rx="10.9091" ry="10.9091" style="fill:#FF7F7F;fill-rule:evenodd;stroke:#FF7F7F;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<ellipse cx="206.364" cy="132.121" rx="10.9091" ry="10.9091" style="fill:#FF7F7F;fill-rule:evenodd;stroke:#FF7F7F;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<ellipse cx="160.909" cy="132.121" rx="10.9091" ry="10.9091" style="fill:#FF7F7F;fill-rule:evenodd;stroke:#FF7F7F;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,147.879 L 206.364,132.121" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,147.879 L 160.909,132.121" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 183.636,147.879 L 183.636,179.394" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 206.364,132.121 L 229.091,147.879" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 160.909,132.121 L 138.182,147.879" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n</svg>',
                                     'width': 400,
                                     'height': 300
                                 }
                             }
                         }
                     }
                 }
             },
             '400': {
                 'description': 'Invalid request data or molecule cannot be processed',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'error': {'type': 'string'}
                             }
                         },
                         'example': {
                             'error': 'Failed to generate molecule visualization'
                         }
                     }
                 }
             },
             '500': {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'error': {'type': 'string'}
                             }
                         },
                         'example': {
                             'error': 'Internal server error'
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(molecule_visualization_fields, code=200, description='Molecule visualization as SVG')
    @marshal_with(molecule_visualization_fields)
    def post(self) -> Tuple[Dict[str, Any], int]:
        """Generate a visualization for a molecule.
        
        Creates an SVG visualization of the provided molecule with customizable
        dimensions and optional atom highlighting. This endpoint is useful for
        visualizing molecular structures in web applications, reports, or
        interactive tools.
        
        Returns:
            Tuple[Dict[str, Any], int]: SVG visualization data and HTTP status code
            
        Raises:
            400: If request data is invalid or molecule cannot be processed
            500: If server error occurs
            
        Request Body:
            molecule_data (str): Molecule data in the specified format (required)
                Example: "C1=CC=C(C=C1)C(=O)O" (Benzoic acid)
            input_format (str): Format of the input data (default: 'smiles')
                Supported formats: 'smiles', 'inchi', 'mol'
            width (int): Width of the SVG in pixels (default: 400)
            height (int): Height of the SVG in pixels (default: 300)
            highlight_atoms (List[int], optional): List of atom indices to highlight (0-based)
                Example: [0, 1, 2] to highlight the first three atoms
            
        Response:
            A JSON object containing:
            - svg (str): SVG string representation of the molecule
            - width (int): Width of the SVG in pixels
            - height (int): Height of the SVG in pixels
            
        Developer Notes:
            - The SVG output can be directly embedded in HTML or used with JavaScript visualization libraries
            - Atom indices are 0-based and correspond to the order in the SMILES/InChI representation
            - For complex molecules, consider increasing the width and height parameters
            - The visualization removes hydrogens for clarity but preserves stereochemistry
            - The SVG can be converted to other formats (PNG, PDF) using standard conversion tools
        """
        try:
            # Validate request data
            schema = MoleculeVisualizationSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Generate SVG
            svg = generate_molecule_svg(
                data["molecule_data"],
                data.get("input_format", "smiles"),
                data.get("width", 400),
                data.get("height", 300),
                data.get("highlight_atoms")
            )
            
            if not svg:
                return _handle_json_serialization({'error': "Failed to generate molecule visualization"}), 400
            
            return _handle_json_serialization({
                "svg": svg,
                "width": data.get("width", 400),
                "height": data.get("height", 300)
            }), 200
        except Exception as e:
            current_app.logger.error(f"Error generating molecule visualization: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500


class SubstructureSearchResource(MethodResource):
    """Resource for performing substructure searches.
    
    This endpoint searches for a specified substructure within a target molecule
    using RDKit's substructure matching capabilities. It can identify specific
    functional groups or structural motifs that may be relevant for
    cryoprotection effectiveness.
    """
    
    @doc(description='Perform a substructure search to find occurrences of a query pattern within a target molecule',
         tags=['RDKit'],
         request_body={
             'query_mol_data': {'description': 'Query substructure data (SMARTS or SMILES pattern to search for)', 'type': 'string', 'required': True, 'example': 'c1ccccc1'},
             'target_mol_data': {'description': 'Target molecule data to search within', 'type': 'string', 'required': True, 'example': 'CC(=O)OC1=CC=CC=C1C(=O)O'},
             'query_format': {'description': 'Format of the query data', 'type': 'string', 'default': 'smarts', 'enum': ['smarts', 'smiles'], 'example': 'smarts'},
             'target_format': {'description': 'Format of the target data', 'type': 'string', 'default': 'smiles', 'enum': ['smiles', 'inchi', 'mol'], 'example': 'smiles'}
         },
         responses={
             '200': {
                 'description': 'Substructure search results',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'match': {'type': 'boolean', 'description': 'Whether a match was found', 'example': True},
                                 'match_count': {'type': 'integer', 'description': 'Number of matches found', 'example': 1},
                                 'matches': {'type': 'array', 'description': 'List of atom indices for each match', 'example': [[0, 1, 2, 3, 4, 5]]},
                                 'visualization': {'type': 'string', 'description': 'SVG visualization with highlighted matches (if found)'}
                             }
                         },
                         'examples': {
                             'Benzene Ring in Aspirin': {
                                 'value': {
                                     'match': True,
                                     'match_count': 1,
                                     'matches': [[1, 2, 3, 4, 5, 6]],
                                     'visualization': '<?xml version="1.0" encoding="UTF-8"?>\n<svg baseProfile="full" height="300px" version="1.1" viewBox="0 0 400 300" width="400px" xml:space="preserve" xmlns="http://www.w3.org/2000/svg" xmlns:rdkit="http://www.rdkit.org/xml" xmlns:xlink="http://www.w3.org/1999/xlink">\n<!-- Truncated SVG content --></svg>'
                                 }
                             },
                             'No Match Found': {
                                 'value': {
                                     'match': False,
                                     'match_count': 0,
                                     'matches': []
                                 }
                             }
                         }
                     }
                 }
             },
             '400': {
                 'description': 'Invalid request data or molecules cannot be processed',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'error': {'type': 'string'}
                             }
                         },
                         'example': {
                             'error': 'Failed to parse molecules'
                         }
                     }
                 }
             },
             '500': {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'error': {'type': 'string'}
                             }
                         },
                         'example': {
                             'error': 'Internal server error'
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(substructure_search_fields, code=200, description='Substructure search results')
    @marshal_with(substructure_search_fields)
    def post(self) -> Tuple[Dict[str, Any], int]:
        """Perform a substructure search.
        
        Searches for a specified substructure (query) within a target molecule
        using RDKit's substructure matching capabilities. Returns match information
        including atom indices and match count.
        
        Returns:
            Tuple[Dict[str, Any], int]: Search results and HTTP status code
            
        Raises:
            400: If request data is invalid or molecules cannot be processed
            500: If server error occurs
            
        Request Body:
            query_mol_data (str): Query substructure data (SMARTS or SMILES pattern) (required)
                Example: "c1ccccc1" (benzene ring)
            target_mol_data (str): Target molecule data to search within (required)
                Example: "CC(=O)OC1=CC=CC=C1C(=O)O" (Aspirin)
            query_format (str): Format of the query data (default: 'smarts')
                Supported formats: 'smarts', 'smiles'
            target_format (str): Format of the target data (default: 'smiles')
                Supported formats: 'smiles', 'inchi', 'mol'
                
        Response:
            Information about matches including:
            - match (bool): Whether a match was found
            - match_count (int): Number of matches found
            - matches (List[List[int]]): List of atom indices for each match
                Example: [[0, 1, 2, 3, 4, 5]] for a single match of a benzene ring
            - visualization (str): SVG visualization with highlighted matches (if found)
            
        Developer Notes:
            - SMARTS is the recommended query format for maximum flexibility in pattern matching
            - For complex substructures, consider using simplified SMARTS patterns
            - The visualization includes highlighted atoms where matches were found
            - Multiple matches will all be highlighted in the same visualization
            - Atom indices are 0-based and correspond to the order in the molecule representation
        """
        try:
            # Validate request data
            schema = SubstructureSearchSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Perform search
            result = perform_substructure_search(
                data["query_mol_data"],
                data["target_mol_data"],
                data.get("query_format", "smarts"),
                data.get("target_format", "smiles")
            )
            
            if "error" in result:
                return _handle_json_serialization({'error': result["error"]}), 400
            
            return _handle_json_serialization(result), 200
        except Exception as e:
            current_app.logger.error(f"Error performing substructure search: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500


class SimilaritySearchResource(MethodResource):
    """Resource for calculating molecular similarity.
    
    This endpoint calculates the similarity between two molecules using various
    fingerprint methods. It helps identify structurally similar compounds that
    may share cryoprotective properties or mechanisms of action.
    """
    
    @doc(description='Calculate similarity between two molecules using various fingerprinting methods',
         tags=['RDKit'],
         request_body={
             'mol1_data': {'description': 'First molecule data', 'type': 'string', 'required': True, 'example': 'CC(=O)OC1=CC=CC=C1C(=O)O'},
             'mol2_data': {'description': 'Second molecule data', 'type': 'string', 'required': True, 'example': 'C1=CC=C(C=C1)C(=O)O'},
             'mol1_format': {'description': 'Format of the first molecule data', 'type': 'string', 'default': 'smiles', 'enum': ['smiles', 'inchi', 'mol'], 'example': 'smiles'},
             'mol2_format': {'description': 'Format of the second molecule data', 'type': 'string', 'default': 'smiles', 'enum': ['smiles', 'inchi', 'mol'], 'example': 'smiles'},
             'fingerprint_type': {'description': 'Type of fingerprint to use for similarity calculation', 'type': 'string', 'default': 'morgan', 'enum': ['morgan', 'maccs', 'rdkit', 'topological'], 'example': 'morgan'}
         },
         responses={
             '200': {
                 'description': 'Molecular similarity results',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'tanimoto': {'type': 'number', 'description': 'Tanimoto similarity coefficient (0-1)', 'example': 0.68},
                                 'dice': {'type': 'number', 'description': 'Dice similarity coefficient (0-1)', 'example': 0.72},
                                 'fingerprint_type': {'type': 'string', 'description': 'Type of fingerprint used for calculation', 'example': 'morgan'}
                             }
                         },
                         'examples': {
                             'Aspirin vs Benzoic Acid': {
                                 'value': {
                                     'tanimoto': 0.68,
                                     'dice': 0.72,
                                     'fingerprint_type': 'morgan'
                                 }
                             },
                             'Similar Compounds': {
                                 'value': {
                                     'tanimoto': 0.92,
                                     'dice': 0.95,
                                     'fingerprint_type': 'morgan'
                                 }
                             }
                         }
                     }
                 }
             },
             '400': {
                 'description': 'Invalid request data or molecules cannot be processed',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'error': {'type': 'string'}
                             }
                         },
                         'example': {
                             'error': 'Failed to parse molecules'
                         }
                     }
                 }
             },
             '500': {
                 'description': 'Server error',
                 'content': {
                     'application/json': {
                         'schema': {
                             'type': 'object',
                             'properties': {
                                 'error': {'type': 'string'}
                             }
                         },
                         'example': {
                             'error': 'Internal server error'
                         }
                     }
                 }
             }
         })
    @apispec_marshal_with(similarity_fields, code=200, description='Molecular similarity results')
    @marshal_with(similarity_fields)
    def post(self) -> Tuple[Dict[str, Any], int]:
        """Calculate similarity between two molecules.
        
        Computes the similarity between two molecules using the specified fingerprint
        method. Returns similarity scores between 0 (completely different) and
        1 (identical).
        
        Returns:
            Tuple[Dict[str, Any], int]: Similarity results and HTTP status code
            
        Raises:
            400: If request data is invalid or molecules cannot be processed
            500: If server error occurs
            
        Request Body:
            mol1_data (str): First molecule data (required)
                Example: "CC(=O)OC1=CC=CC=C1C(=O)O" (Aspirin)
            mol2_data (str): Second molecule data (required)
                Example: "C1=CC=C(C=C1)C(=O)O" (Benzoic acid)
            mol1_format (str): Format of the first molecule data (default: 'smiles')
                Supported formats: 'smiles', 'inchi', 'mol'
            mol2_format (str): Format of the second molecule data (default: 'smiles')
                Supported formats: 'smiles', 'inchi', 'mol'
            fingerprint_type (str): Type of fingerprint to use (default: 'morgan')
                Supported types: 'morgan', 'maccs', 'rdkit', 'topological'
                
        Response:
            Similarity information including:
            - tanimoto (float): Tanimoto similarity coefficient (0-1)
                Example: 0.68 for moderately similar molecules
            - dice (float): Dice similarity coefficient (0-1)
                Example: 0.72 for moderately similar molecules
            - fingerprint_type (str): Type of fingerprint used for calculation
                Example: "morgan"
                
        Developer Notes:
            - Morgan (ECFP) fingerprints are recommended for most similarity calculations
            - MACCS keys are useful for pharmacophore-based similarity
            - Tanimoto coefficient is the most widely used similarity metric
            - Dice coefficient gives more weight to common features than Tanimoto
            - Similarity values above 0.85 typically indicate highly similar compounds
            - Values below 0.4 typically indicate structurally distinct compounds
        """
        try:
            # Validate request data
            schema = SimilaritySearchSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                return _handle_json_serialization({'error': str(err.messages)}), 400
            
            # Calculate similarity
            result = calculate_similarity(
                data["mol1_data"],
                data["mol2_data"],
                data.get("mol1_format", "smiles"),
                data.get("mol2_format", "smiles"),
                data.get("fingerprint_type", "morgan")
            )
            
            if "error" in result:
                return _handle_json_serialization({'error': result["error"]}), 400
            
            return _handle_json_serialization(result), 200
        except Exception as e:
            current_app.logger.error(f"Error calculating molecular similarity: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500


class MoleculePropertyCalculationResource(MethodResource):
    """Resource for calculating properties for a specific molecule in the database.
    
    This endpoint calculates properties for a molecule that already exists in the
    database and stores the results. It provides a convenient way to enrich the
    database with computed molecular properties for further analysis.
    """
    
    @doc(description='Calculate and store properties for a molecule in the database',
         tags=['RDKit'],
         params={'molecule_id': {'description': 'ID of the molecule to calculate properties for', 'type': 'string', 'required': True}},
         security=[{'Bearer': []}])
    @token_required
    def post(self, molecule_id: str) -> Tuple[Dict[str, Any], int]:
        """Calculate and store properties for a molecule in the database.
        
        Retrieves a molecule from the database by ID, calculates its properties
        using RDKit, and stores the results in the molecular_properties table.
        
        Args:
            molecule_id (str): ID of the molecule to calculate properties for
            
        Returns:
            Tuple[Dict[str, Any], int]: Success message and HTTP status code
            
        Raises:
            400: If molecule data is invalid or properties cannot be calculated
            404: If molecule not found
            500: If server error occurs
            
        Authentication:
            Requires a valid JWT token in the Authorization header.
            Example: Authorization: Bearer <token>
            
        Response:
            Success message with the number of properties calculated and stored
        """
        try:
            from api.utils import get_supabase_client, get_user_id
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Get the molecule from the database
            response = supabase.table("molecules").select("*").eq("id", molecule_id).execute()
            
            if not response.data:
                return _handle_json_serialization({'error': f"Molecule with ID {molecule_id} not found"}), 404
            
            molecule = response.data[0]
            smiles = molecule.get("smiles")
            
            if not smiles:
                return _handle_json_serialization({'error': "Molecule does not have SMILES data"}), 400
            
            # Calculate properties
            properties = calculate_all_properties(smiles, "smiles")
            
            if "error" in properties:
                return _handle_json_serialization({'error': properties["error"]}), 400
            
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
                    return _handle_json_serialization({'error': f"Error inserting properties: {response.error}"}), 400
                
                return _handle_json_serialization({"message": f"Calculated and stored {len(property_inserts)} properties for molecule {molecule_id}"}), 200
            else:
                return _handle_json_serialization({"message": "No properties to store"}), 200
            
        except Exception as e:
            current_app.logger.error(f"Error calculating properties for molecule: {str(e)}")
            return _handle_json_serialization(handle_supabase_error(e)), 500