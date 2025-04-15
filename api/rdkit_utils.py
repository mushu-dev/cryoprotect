"""
CryoProtect Analyzer API - RDKit Utilities

This module provides functions for molecular property calculations using RDKit.
It includes functions for calculating hydrogen bonding capacity, XLogP,
topological polar surface area, molecular weight and volume, functional group
identification, and permeability coefficient estimation.
"""

import logging
from typing import Dict, List, Optional, Tuple, Union, Any

# Import RDKit modules
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, MolSurf, AllChem, Fragments
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    
    # Try to import IPythonConsole, but don't fail if it's not available
    try:
        from rdkit.Chem.Draw import IPythonConsole
        IPYTHON_AVAILABLE = True
    except ImportError:
        logging.warning("IPython is not installed. Some visualization features may be limited.")
        IPYTHON_AVAILABLE = False
except ImportError:
    logging.error("RDKit is not installed. Please install it using conda: conda install -c conda-forge rdkit")
    raise

# Set up logging
logger = logging.getLogger(__name__)

def parse_molecule(mol_data: str, input_format: str = 'smiles') -> Optional[Chem.Mol]:
    """
    Parse molecular data in various formats into an RDKit molecule object.
    
    Args:
        mol_data: Molecular data as a string (SMILES, MOL, SDF)
        input_format: Format of the input data ('smiles', 'mol', 'sdf')
        
    Returns:
        RDKit Mol object or None if parsing fails
    """
    try:
        if input_format.lower() == 'smiles':
            mol = Chem.MolFromSmiles(mol_data)
        elif input_format.lower() == 'mol':
            mol = Chem.MolFromMolBlock(mol_data)
        elif input_format.lower() == 'sdf':
            # For SDF, we take the first molecule
            suppl = Chem.SDMolSupplier()
            suppl.SetData(mol_data)
            mol = next(suppl) if suppl else None
        else:
            logger.error(f"Unsupported input format: {input_format}")
            return None
            
        if mol is None:
            logger.error(f"Failed to parse molecule from {input_format} data")
            return None
            
        # Add hydrogens and compute 3D coordinates if needed
        mol = Chem.AddHs(mol)
        if not mol.GetNumConformers():
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
        return mol
    except Exception as e:
        logger.error(f"Error parsing molecule: {str(e)}")
        return None

def calculate_hydrogen_bonding(mol: Chem.Mol) -> Dict[str, int]:
    """
    Calculate hydrogen bond donors and acceptors.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        Dictionary with donor and acceptor counts
    """
    if mol is None:
        return {"donors": 0, "acceptors": 0}
        
    donors = Lipinski.NumHDonors(mol)
    acceptors = Lipinski.NumHAcceptors(mol)
    
    return {
        "donors": donors,
        "acceptors": acceptors,
        "total": donors + acceptors
    }

def calculate_logp(mol: Chem.Mol) -> float:
    """
    Calculate XLogP (partition coefficient).
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        XLogP value
    """
    if mol is None:
        return 0.0
        
    return Descriptors.MolLogP(mol)

def calculate_tpsa(mol: Chem.Mol) -> float:
    """
    Calculate topological polar surface area.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        TPSA value in Å²
    """
    if mol is None:
        return 0.0
        
    return Descriptors.TPSA(mol)

def calculate_molecular_properties(mol: Chem.Mol) -> Dict[str, float]:
    """
    Calculate basic molecular properties.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        Dictionary of molecular properties
    """
    if mol is None:
        return {}
        
    properties = {
        "molecular_weight": Descriptors.MolWt(mol),
        "exact_mass": Descriptors.ExactMolWt(mol),
        "heavy_atom_count": Descriptors.HeavyAtomCount(mol),
        "atom_count": mol.GetNumAtoms(),
        "rotatable_bond_count": Descriptors.NumRotatableBonds(mol),
        "ring_count": Descriptors.RingCount(mol),
        "aromatic_ring_count": Lipinski.NumAromaticRings(mol),
        "fraction_csp3": Descriptors.FractionCSP3(mol)
    }
    
    # Calculate molecular volume (requires 3D coordinates)
    if mol.GetNumConformers() > 0:
        properties["molecular_volume"] = AllChem.ComputeMolVolume(mol)
    
    return properties

def identify_functional_groups(mol: Chem.Mol) -> Dict[str, int]:
    """
    Identify and count functional groups in the molecule.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        Dictionary of functional groups and their counts
    """
    if mol is None:
        return {}
        
    functional_groups = {
        # Alcohols (using available fragments)
        "alcohol": Fragments.fr_Al_OH(mol) + Fragments.fr_Ar_OH(mol),
        "aldehyde": Fragments.fr_aldehyde(mol),
        "alkyl_halide": Fragments.fr_alkyl_halide(mol),
        "amide": Fragments.fr_amide(mol),
        # Use NH groups instead of generic amine
        "amine": Fragments.fr_NH1(mol) + Fragments.fr_NH2(mol),
        "carboxylic_acid": Fragments.fr_COO(mol),
        "ester": Fragments.fr_ester(mol),
        "ether": Fragments.fr_ether(mol),
        "ketone": Fragments.fr_ketone(mol),
        "nitrile": Fragments.fr_nitrile(mol),
        "nitro": Fragments.fr_nitro(mol),
        "sulfide": Fragments.fr_sulfide(mol),
        # Use sulfonamd instead of sulfonamide
        "sulfonamide": Fragments.fr_sulfonamd(mol),
        "sulfone": Fragments.fr_sulfone(mol),
        # Remove sulfonic_acid as fr_sulfonic doesn't exist
        "phenol": Fragments.fr_phenol(mol),
        # Remove phosphate and phosphonate as they don't exist
        # Use combined hydroxyl groups
        "hydroxyl": Fragments.fr_Al_OH(mol) + Fragments.fr_Ar_OH(mol)
    }
    
    # Filter out groups with zero count
    return {k: v for k, v in functional_groups.items() if v > 0}

def estimate_permeability(mol: Chem.Mol) -> Dict[str, float]:
    """
    Estimate permeability coefficients based on molecular properties.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        Dictionary of permeability estimates
    """
    if mol is None:
        return {}
        
    # Calculate properties relevant to permeability
    logp = calculate_logp(mol)
    tpsa = calculate_tpsa(mol)
    mw = Descriptors.MolWt(mol)
    hbond_donors = Lipinski.NumHDonors(mol)
    hbond_acceptors = Lipinski.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    # Estimate cell membrane permeability using various models
    
    # Rule of 5 violations (Lipinski's rule)
    rule_of_5_violations = sum([
        logp > 5,
        mw > 500,
        hbond_donors > 5,
        hbond_acceptors > 10
    ])
    
    # Veber rule violations
    veber_violations = sum([
        rotatable_bonds > 10,
        tpsa > 140
    ])
    
    # BOILED-Egg model for blood-brain barrier (BBB) and intestinal absorption
    # Reference: Daina & Zoete, 2016
    bbb_permeant = (logp < 3 and tpsa < 75) or (logp > 3 and tpsa < 60)
    intestinal_absorption = (logp < 5 and tpsa < 130)
    
    # Estimate permeability coefficient (log Papp) using a simple model
    # This is a simplified model based on TPSA and LogP
    log_papp = -4.3 - 0.01 * tpsa + 0.4 * logp
    
    return {
        "rule_of_5_violations": rule_of_5_violations,
        "veber_violations": veber_violations,
        "bbb_permeant": bbb_permeant,
        "intestinal_absorption": intestinal_absorption,
        "estimated_log_papp": log_papp
    }

def calculate_all_properties(mol_data: str, input_format: str = 'smiles') -> Dict[str, Any]:
    """
    Calculate all molecular properties for a given molecule.
    
    Args:
        mol_data: Molecular data as a string (SMILES, MOL, SDF)
        input_format: Format of the input data ('smiles', 'mol', 'sdf')
        
    Returns:
        Dictionary of all calculated properties
    """
    mol = parse_molecule(mol_data, input_format)
    
    if mol is None:
        return {"error": "Failed to parse molecule"}
    
    # Calculate all properties
    properties = {}
    properties.update({"hydrogen_bonding": calculate_hydrogen_bonding(mol)})
    properties.update({"logp": calculate_logp(mol)})
    properties.update({"tpsa": calculate_tpsa(mol)})
    properties.update({"molecular_properties": calculate_molecular_properties(mol)})
    properties.update({"functional_groups": identify_functional_groups(mol)})
    properties.update({"permeability": estimate_permeability(mol)})
    
    # Add SMILES and InChI for reference
    properties["smiles"] = Chem.MolToSmiles(mol)
    properties["inchi"] = Chem.MolToInchi(mol)
    properties["inchi_key"] = Chem.MolToInchiKey(mol)
    
    return properties

def generate_molecule_svg(mol_data: str, input_format: str = 'smiles', 
                         width: int = 400, height: int = 300, 
                         highlight_atoms: List[int] = None) -> str:
    """
    Generate an SVG image of a molecule.
    
    Args:
        mol_data: Molecular data as a string (SMILES, MOL, SDF)
        input_format: Format of the input data ('smiles', 'mol', 'sdf')
        width: Image width in pixels
        height: Image height in pixels
        highlight_atoms: List of atom indices to highlight
        
    Returns:
        SVG string representation of the molecule
    """
    mol = parse_molecule(mol_data, input_format)
    
    if mol is None:
        return ""
    
    # Remove hydrogens for cleaner visualization
    mol = Chem.RemoveHs(mol)
    
    # Create a drawing object
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    
    # Set drawing options
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.addAtomIndices = False
    
    # Draw the molecule
    if highlight_atoms:
        drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
    else:
        drawer.DrawMolecule(mol)
    
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    
    return svg

def perform_substructure_search(query_mol_data: str, target_mol_data: str, 
                               query_format: str = 'smarts', target_format: str = 'smiles') -> Dict[str, Any]:
    """
    Perform a substructure search.
    
    Args:
        query_mol_data: Query molecule data (usually SMARTS pattern)
        target_mol_data: Target molecule data to search in
        query_format: Format of the query data ('smarts', 'smiles')
        target_format: Format of the target data ('smiles', 'mol', 'sdf')
        
    Returns:
        Dictionary with search results
    """
    # Parse query molecule
    if query_format.lower() == 'smarts':
        query_mol = Chem.MolFromSmarts(query_mol_data)
    else:
        query_mol = parse_molecule(query_mol_data, query_format)
    
    # Parse target molecule
    target_mol = parse_molecule(target_mol_data, target_format)
    
    if query_mol is None or target_mol is None:
        return {"match": False, "error": "Failed to parse molecules"}
    
    # Perform substructure search
    matches = target_mol.GetSubstructMatches(query_mol)
    
    result = {
        "match": len(matches) > 0,
        "match_count": len(matches),
        "matches": [list(match) for match in matches]
    }
    
    # Generate visualization if there are matches
    if matches:
        # Generate SVG with highlighted matches
        all_match_atoms = set()
        for match in matches:
            all_match_atoms.update(match)
        
        result["visualization"] = generate_molecule_svg(
            Chem.MolToSmiles(target_mol), 
            'smiles', 
            highlight_atoms=list(all_match_atoms)
        )
    
    return result

def calculate_similarity(mol1_data: str, mol2_data: str, 
                        mol1_format: str = 'smiles', mol2_format: str = 'smiles',
                        fingerprint_type: str = 'morgan') -> Dict[str, float]:
    """
    Calculate molecular similarity between two molecules.
    
    Args:
        mol1_data: First molecule data
        mol2_data: Second molecule data
        mol1_format: Format of the first molecule data
        mol2_format: Format of the second molecule data
        fingerprint_type: Type of fingerprint to use ('morgan', 'maccs', 'topological')
        
    Returns:
        Dictionary with similarity metrics
    """
    from rdkit import DataStructs
    from rdkit.Chem import AllChem, MACCSkeys
    
    # Parse molecules
    mol1 = parse_molecule(mol1_data, mol1_format)
    mol2 = parse_molecule(mol2_data, mol2_format)
    
    if mol1 is None or mol2 is None:
        return {"error": "Failed to parse molecules"}
    
    # Generate fingerprints based on the specified type
    if fingerprint_type.lower() == 'morgan':
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    elif fingerprint_type.lower() == 'maccs':
        fp1 = MACCSkeys.GenMACCSKeys(mol1)
        fp2 = MACCSkeys.GenMACCSKeys(mol2)
    elif fingerprint_type.lower() == 'topological':
        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)
    else:
        return {"error": f"Unsupported fingerprint type: {fingerprint_type}"}
    
    # Calculate similarity metrics
    tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)
    dice = DataStructs.DiceSimilarity(fp1, fp2)
    
    return {
        "tanimoto": tanimoto,
        "dice": dice,
        "fingerprint_type": fingerprint_type
    }