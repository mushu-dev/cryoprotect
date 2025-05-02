"""
CryoProtect Analyzer API - RDKit Schemas

This module contains Marshmallow schemas for RDKit-related API endpoints.
These schemas are used for request validation and documentation generation.
"""

from marshmallow import Schema, fields, validate

class MoleculeInputSchema(Schema):
    """
    Schema for molecule input data.
    
    This schema validates the input data for a molecule, including the molecular
    structure data and its format. It ensures that all required fields are present
    and that values meet validation requirements.
    
    Attributes:
        molecule_data: Molecular structure data (required)
        input_format: Format of the molecule data (default: 'smiles')
    
    Validation:
        - Input format must be one of: smiles, inchi, mol
    """
    molecule_data = fields.String(
        required=True,
        description="Molecular structure data (SMILES, InChI, or MOL format)",
        example="CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin SMILES
    )
    input_format = fields.String(
        validate=validate.OneOf(['smiles', 'inchi', 'mol']),
        default='smiles',
        description="Format of the molecule data",
        example="smiles"
    )

class MoleculeVisualizationSchema(Schema):
    """
    Schema for molecule visualization options.
    
    This schema validates the options for visualizing a molecule, including the
    molecular structure data, format, dimensions, and atom highlighting. It ensures
    that all required fields are present and that values meet validation requirements.
    
    Attributes:
        molecule_data: Molecular structure data (required)
        input_format: Format of the molecule data (default: 'smiles')
        width: Width of the visualization in pixels (default: 400)
        height: Height of the visualization in pixels (default: 300)
        highlight_atoms: List of atom indices to highlight (optional)
    
    Validation:
        - Input format must be one of: smiles, inchi, mol
        - Width and height must be positive integers
    """
    molecule_data = fields.String(
        required=True,
        description="Molecular structure data (SMILES, InChI, or MOL format)",
        example="C1=CC=C(C=C1)C(=O)O"  # Benzoic acid SMILES
    )
    input_format = fields.String(
        validate=validate.OneOf(['smiles', 'inchi', 'mol']),
        default='smiles',
        description="Format of the molecule data",
        example="smiles"
    )
    width = fields.Integer(
        default=400,
        validate=validate.Range(min=1),
        description="Width of the visualization in pixels",
        example=400
    )
    height = fields.Integer(
        default=300,
        validate=validate.Range(min=1),
        description="Height of the visualization in pixels",
        example=300
    )
    highlight_atoms = fields.List(
        fields.Integer(),
        description="List of atom indices to highlight (0-based)",
        example=[0, 1, 2]
    )

class SubstructureSearchSchema(Schema):
    """
    Schema for substructure search.
    
    This schema validates the parameters for a substructure search, including
    the query and target molecule data and their formats. It ensures that all
    required fields are present and that values meet validation requirements.
    
    Attributes:
        query_mol_data: Query molecule data (SMARTS or SMILES pattern) (required)
        target_mol_data: Target molecule data to search within (required)
        query_format: Format of the query molecule data (default: 'smarts')
        target_format: Format of the target molecule data (default: 'smiles')
    
    Validation:
        - Query format must be one of: smarts, smiles
        - Target format must be one of: smiles, inchi, mol
    """
    query_mol_data = fields.String(
        required=True,
        description="Query molecule data (SMARTS or SMILES pattern to search for)",
        example="c1ccccc1"  # Benzene ring SMARTS
    )
    target_mol_data = fields.String(
        required=True,
        description="Target molecule data to search within",
        example="CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin SMILES
    )
    query_format = fields.String(
        validate=validate.OneOf(['smarts', 'smiles']),
        default='smarts',
        description="Format of the query molecule data",
        example="smarts"
    )
    target_format = fields.String(
        validate=validate.OneOf(['smiles', 'inchi', 'mol']),
        default='smiles',
        description="Format of the target molecule data",
        example="smiles"
    )

class SimilaritySearchSchema(Schema):
    """
    Schema for similarity search.
    
    This schema validates the parameters for a molecular similarity search,
    including the two molecules to compare, their formats, and the fingerprint
    type to use for the comparison. It ensures that all required fields are present
    and that values meet validation requirements.
    
    Attributes:
        mol1_data: First molecule data (required)
        mol2_data: Second molecule data to compare with (required)
        mol1_format: Format of the first molecule data (default: 'smiles')
        mol2_format: Format of the second molecule data (default: 'smiles')
        fingerprint_type: Type of molecular fingerprint to use (default: 'morgan')
    
    Validation:
        - Molecule formats must be one of: smiles, inchi, mol
        - Fingerprint type must be one of: morgan, maccs, rdkit, topological
    """
    mol1_data = fields.String(
        required=True,
        description="First molecule data",
        example="CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin SMILES
    )
    mol2_data = fields.String(
        required=True,
        description="Second molecule data to compare with",
        example="CC(=O)OC1=CC=CC=C1"  # Acetylphenol SMILES
    )
    mol1_format = fields.String(
        validate=validate.OneOf(['smiles', 'inchi', 'mol']),
        default='smiles',
        description="Format of the first molecule data",
        example="smiles"
    )
    mol2_format = fields.String(
        validate=validate.OneOf(['smiles', 'inchi', 'mol']),
        default='smiles',
        description="Format of the second molecule data",
        example="smiles"
    )
    fingerprint_type = fields.String(
        validate=validate.OneOf(['morgan', 'maccs', 'rdkit', 'topological']),
        default='morgan',
        description="Type of molecular fingerprint to use for similarity calculation",
        example="morgan"
    )

class BatchMoleculeSchema(Schema):
    """
    Schema for batch molecule processing.
    
    This schema validates the data for batch molecule processing, including a list
    of molecules with their data and formats. It ensures that all required fields
    are present and that values meet validation requirements.
    
    Attributes:
        molecules: List of molecule data objects (required)
    
    Validation:
        - At least one molecule must be provided
        - Each molecule must have data and format fields
    """
    molecules = fields.List(
        fields.Dict(keys=fields.String(), values=fields.String()),
        required=True,
        validate=validate.Length(min=1),
        description="List of molecule data objects, each with 'data' and 'format' fields",
        example=[
            {"data": "CC(=O)OC1=CC=CC=C1C(=O)O", "format": "smiles"},
            {"data": "C1=CC=C(C=C1)C(=O)O", "format": "smiles"}
        ]
    )

class SimilaritySearchBatchSchema(Schema):
    """
    Schema for batch similarity search.
    
    This schema validates the parameters for a batch similarity search, including
    the query molecule and a list of target molecules to compare against. It ensures
    that all required fields are present and that values meet validation requirements.
    
    Attributes:
        query_mol_data: Query molecule data (required)
        target_molecules: List of target molecule data strings (required)
        query_format: Format of the query molecule data (default: 'smiles')
        fingerprint_type: Type of molecular fingerprint to use (default: 'morgan')
        similarity_threshold: Minimum similarity score to include in results (default: 0.7)
    
    Validation:
        - Query format must be one of: smiles, mol, sdf
        - Fingerprint type must be one of: morgan, maccs, topological
        - Similarity threshold must be between 0 and 1
    """
    query_mol_data = fields.String(
        required=True,
        description="Query molecule data",
        example="CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin SMILES
    )
    target_molecules = fields.List(
        fields.String(),
        required=True,
        validate=validate.Length(min=1),
        description="List of target molecule data strings to compare against",
        example=["C1=CC=C(C=C1)C(=O)O", "CC(=O)OC1=CC=CC=C1"]
    )
    query_format = fields.String(
        validate=validate.OneOf(['smiles', 'mol', 'sdf']),
        default='smiles',
        description="Format of the query molecule data",
        example="smiles"
    )
    fingerprint_type = fields.String(
        validate=validate.OneOf(['morgan', 'maccs', 'topological']),
        default='morgan',
        description="Type of molecular fingerprint to use for similarity calculation",
        example="morgan"
    )
    similarity_threshold = fields.Float(
        validate=validate.Range(min=0, max=1),
        default=0.7,
        description="Minimum similarity score to include in results (0-1)",
        example=0.7
    )

class SubstructureSearchBatchSchema(Schema):
    """
    Schema for batch substructure search.
    
    This schema validates the parameters for a batch substructure search, including
    the query substructure and a list of target molecules to search within. It ensures
    that all required fields are present and that values meet validation requirements.
    
    Attributes:
        query_mol_data: Query substructure data (SMARTS or SMILES pattern) (required)
        target_molecules: List of target molecule data strings (required)
        query_format: Format of the query substructure data (default: 'smarts')
    
    Validation:
        - Query format must be one of: smarts, smiles
    """
    query_mol_data = fields.String(
        required=True,
        description="Query substructure data (SMARTS or SMILES pattern to search for)",
        example="c1ccccc1"  # Benzene ring SMARTS
    )
    target_molecules = fields.List(
        fields.String(),
        required=True,
        validate=validate.Length(min=1),
        description="List of target molecule data strings to search within",
        example=["CC(=O)OC1=CC=CC=C1C(=O)O", "C1=CC=C(C=C1)C(=O)O"]
    )
    query_format = fields.String(
        validate=validate.OneOf(['smarts', 'smiles']),
        default='smarts',
        description="Format of the query substructure data",
        example="smarts"
    )

class MoleculeGridSchema(Schema):
    """
    Schema for molecule grid visualization.
    
    This schema validates the parameters for generating a grid of molecule visualizations,
    including the list of molecules, optional labels, and layout options. It ensures
    that all required fields are present and that values meet validation requirements.
    
    Attributes:
        molecules: List of molecule data strings (required)
        labels: Optional list of labels for each molecule
        mol_per_row: Number of molecules per row in the grid (default: 3)
        sub_img_width: Width of each molecule image in pixels (default: 200)
        sub_img_height: Height of each molecule image in pixels (default: 200)
        input_format: Format of the molecule data (default: 'smiles')
    
    Validation:
        - Input format must be one of: smiles, mol, sdf
        - At least one molecule must be provided
    """
    molecules = fields.List(
        fields.String(),
        required=True,
        validate=validate.Length(min=1),
        description="List of molecule data strings",
        example=["CC(=O)OC1=CC=CC=C1C(=O)O", "C1=CC=C(C=C1)C(=O)O"]
    )
    labels = fields.List(
        fields.String(),
        description="Optional list of labels for each molecule",
        example=["Aspirin", "Benzoic acid"]
    )
    mol_per_row = fields.Integer(
        default=3,
        validate=validate.Range(min=1),
        description="Number of molecules per row in the grid",
        example=3
    )
    sub_img_width = fields.Integer(
        default=200,
        validate=validate.Range(min=1),
        description="Width of each molecule image in pixels",
        example=200
    )
    sub_img_height = fields.Integer(
        default=200,
        validate=validate.Range(min=1),
        description="Height of each molecule image in pixels",
        example=200
    )
    input_format = fields.String(
        validate=validate.OneOf(['smiles', 'mol', 'sdf']),
        default='smiles',
        description="Format of the molecule data",
        example="smiles"
    )

class ScaffoldAnalysisSchema(Schema):
    """
    Schema for scaffold analysis.
    
    This schema validates the parameters for analyzing the scaffold of a molecule,
    including the molecule data and its format. It ensures that all required fields
    are present and that values meet validation requirements.
    
    Attributes:
        molecule_data: Molecular structure data (required)
        input_format: Format of the molecule data (default: 'smiles')
    
    Validation:
        - Input format must be one of: smiles, mol, sdf
    """
    molecule_data = fields.String(
        required=True,
        description="Molecular structure data (SMILES, MOL, or SDF format)",
        example="CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin SMILES
    )
    input_format = fields.String(
        validate=validate.OneOf(['smiles', 'mol', 'sdf']),
        default='smiles',
        description="Format of the molecule data",
        example="smiles"
    )

class MolecularDynamicsSchema(Schema):
    """
    Schema for molecular dynamics simulation.
    
    This schema validates the parameters for a molecular dynamics simulation,
    including the molecule data, simulation parameters, and output options. It ensures
    that all required fields are present and that values meet validation requirements.
    
    Attributes:
        molecule_data: Molecular structure data (required)
        input_format: Format of the molecule data (default: 'smiles')
        simulation_time: Simulation time in picoseconds (default: 10.0)
        temperature: Simulation temperature in Kelvin (default: 300.0)
        force_field: Force field to use for the simulation (default: 'MMFF94')
        solvent: Solvent environment for the simulation (default: 'water')
        include_trajectory: Whether to include trajectory data in the results (default: False)
    
    Validation:
        - Input format must be one of: smiles, mol, sdf
        - Force field must be one of: MMFF94, UFF, AMBER
    """
    molecule_data = fields.String(
        required=True,
        description="Molecular structure data (SMILES, MOL, or SDF format)",
        example="CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin SMILES
    )
    input_format = fields.String(
        validate=validate.OneOf(['smiles', 'mol', 'sdf']),
        default='smiles',
        description="Format of the molecule data",
        example="smiles"
    )
    simulation_time = fields.Float(
        default=10.0,
        validate=validate.Range(min=0.1),
        description="Simulation time in picoseconds",
        example=10.0
    )
    temperature = fields.Float(
        default=300.0,
        validate=validate.Range(min=0.1),
        description="Simulation temperature in Kelvin",
        example=300.0
    )
    force_field = fields.String(
        validate=validate.OneOf(['MMFF94', 'UFF', 'AMBER']),
        default='MMFF94',
        description="Force field to use for the simulation",
        example="MMFF94"
    )
    solvent = fields.String(
        default='water',
        description="Solvent environment for the simulation",
        example="water"
    )
    include_trajectory = fields.Boolean(
        default=False,
        description="Whether to include trajectory data in the results",
        example=False
    )

class MoleculePropertyCalculationBatchSchema(Schema):
    """
    Schema for batch calculation of properties for molecules in the database.
    
    This schema validates the parameters for calculating properties for multiple
    molecules in the database. It ensures that the required fields are present
    and that values meet validation requirements.
    
    Attributes:
        molecule_ids: List of molecule IDs to calculate properties for (required)
    
    Validation:
        - At least one molecule ID must be provided
    """
    molecule_ids = fields.List(
        fields.UUID(),
        required=True,
        validate=validate.Length(min=1),
        description="List of molecule IDs to calculate properties for",
        example=["123e4567-e89b-12d3-a456-426614174000", "223e4567-e89b-12d3-a456-426614174001"]
    )