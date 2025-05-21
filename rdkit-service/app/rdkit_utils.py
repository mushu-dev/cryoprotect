"""
CryoProtect Analyzer API - RDKit Utilities

This module provides core scientific utilities for molecular property calculations using RDKit.

Overview:
    - Implements algorithms for calculating hydrogen bonding capacity, XLogP (partition coefficient),
      topological polar surface area (TPSA), molecular weight, molecular volume, and other key descriptors.
    - Identifies and counts functional groups using RDKit fragment analysis.
    - Estimates permeability coefficients using established medicinal chemistry rules (e.g., Lipinski's Rule of 5, Veber rule)
      and the BOILED-Egg model for blood-brain barrier and intestinal absorption prediction.
    - Provides substructure search and molecular similarity calculations using various fingerprinting methods.

Scientific Context:
    These calculations are foundational in cheminformatics and drug design, as they relate to molecular bioavailability,
    permeability, and pharmacokinetic properties. The algorithms implemented here are based on widely accepted
    literature, including:
        - Lipinski, C. A., et al. (1997). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Advanced Drug Delivery Reviews, 23(1-3), 3-25.
        - Veber, D. F., et al. (2002). Molecular properties that influence the oral bioavailability of drug candidates. Journal of Medicinal Chemistry, 45(12), 2615-2623.
        - Daina, A., & Zoete, V. (2016). A BOILED-Egg to predict gastrointestinal absorption and brain penetration of small molecules. ChemMedChem, 11(11), 1117-1121.

All functions and algorithms are designed to be accessible and maintainable for both developers and scientists.
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
        # Use module-level logger after it's defined
        IPYTHON_AVAILABLE = False
except ImportError:
    # Use module-level logger after it's defined
    raise

# Set up logging
logger = logging.getLogger(__name__)

def parse_molecule(mol_data: str, input_format: str = 'smiles') -> Optional[Chem.Mol]:
    """
    Parses molecular data in various formats into an RDKit molecule object.

    Scientific Rationale:
        Molecular parsing is the first step in cheminformatics workflows.
        SMILES, MOL, and SDF are standard formats for representing chemical structures.
        Hydrogens are added and 3D coordinates are generated to enable property calculations that require 3D geometry.

    Args:
        mol_data (str): Molecular data as a string (SMILES, MOL, SDF).
        input_format (str): Format of the input data ('smiles', 'mol', 'sdf').

    Returns:
        Optional[Chem.Mol]: RDKit Mol object or None if parsing fails.

    References:
        - Weininger, D. (1988). SMILES, a chemical language and information system. 1. Introduction to methodology and encoding rules. J. Chem. Inf. Comput. Sci., 28(1), 31-36.
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
    Calculates the number of hydrogen bond donors and acceptors in a molecule.

    Scientific Rationale:
        Hydrogen bonding capacity is a key determinant of solubility, permeability, and bioavailability.
        Lipinski's Rule of 5 uses these counts to predict drug-likeness.

    Args:
        mol (Chem.Mol): RDKit Mol object.

    Returns:
        Dict[str, int]: Dictionary with donor, acceptor, and total counts.

    References:
        - Lipinski, C. A., et al. (1997). Advanced Drug Delivery Reviews, 23(1-3), 3-25.
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
    Calculates the XLogP (octanol-water partition coefficient) of a molecule.

    Scientific Rationale:
        LogP is a measure of hydrophobicity and is critical for predicting membrane permeability and solubility.
        High logP values indicate lipophilicity, which can affect absorption and distribution.

    Args:
        mol (Chem.Mol): RDKit Mol object.

    Returns:
        float: XLogP value.

    References:
        - Mannhold, R., et al. (2009). Calculation of molecular lipophilicity: State-of-the-art and comparison of logP methods on more than 96,000 compounds. J. Pharm. Sci., 98(3), 861-893.
    """
    if mol is None:
        return 0.0
        
    return Descriptors.MolLogP(mol)

def calculate_tpsa(mol: Chem.Mol) -> float:
    """
    Calculates the topological polar surface area (TPSA) of a molecule.

    Scientific Rationale:
        TPSA is a predictor of drug transport properties, including intestinal absorption and blood-brain barrier penetration.
        Molecules with high TPSA tend to have lower membrane permeability.

    Args:
        mol (Chem.Mol): RDKit Mol object.

    Returns:
        float: TPSA value in Å².

    References:
        - Ertl, P., et al. (2000). Fast calculation of molecular polar surface area as a sum of fragment-based contributions and its application to the prediction of drug transport properties. J. Med. Chem., 43(20), 3714-3717.
    """
    if mol is None:
        return 0.0
        
    return Descriptors.TPSA(mol)

def calculate_molecular_properties(mol: Chem.Mol) -> Dict[str, float]:
    """
    Calculates a set of fundamental molecular properties.

    Scientific Rationale:
        These properties (molecular weight, heavy atom count, rotatable bonds, etc.) are used in medicinal chemistry
        to assess drug-likeness, synthetic accessibility, and physical behavior.

    Args:
        mol (Chem.Mol): RDKit Mol object.

    Returns:
        Dict[str, float]: Dictionary of molecular properties.

    Notes:
        - Molecular volume is calculated only if 3D coordinates are available.
        - FractionCSP3 is a measure of saturation, relevant for oral bioavailability.

    References:
        - Lovering, F., et al. (2009). Escape from flatland: increasing saturation as an approach to improving clinical success. J. Med. Chem., 52(21), 6752-6756.
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
    Identifies and counts key functional groups in the molecule using RDKit fragment analysis.

    Scientific Rationale:
        Functional groups determine chemical reactivity, solubility, and biological activity.
        This function uses RDKit's fragment-based approach to identify groups relevant to cryoprotectant and drug design.

    Args:
        mol (Chem.Mol): RDKit Mol object.

    Returns:
        Dict[str, int]: Dictionary of functional groups and their counts.

    Notes:
        - Only groups with nonzero counts are returned.
        - The selection of groups is based on their relevance to cryoprotectant chemistry.

    References:
        - RDKit documentation: https://www.rdkit.org/docs/source/rdkit.Chem.Fragments.html
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
    Estimates permeability coefficients and related properties using established medicinal chemistry rules.

    Scientific Rationale:
        - Lipinski's Rule of 5 and Veber rule are used to predict oral bioavailability and permeability.
        - The BOILED-Egg model predicts blood-brain barrier (BBB) and intestinal absorption based on logP and TPSA.
        - The estimated log Papp (apparent permeability) is a simplified model for passive membrane diffusion.

    Args:
        mol (Chem.Mol): RDKit Mol object.

    Returns:
        Dict[str, float]: Dictionary of permeability estimates and rule violations.

    References:
        - Lipinski, C. A., et al. (1997). Advanced Drug Delivery Reviews, 23(1-3), 3-25.
        - Veber, D. F., et al. (2002). J. Med. Chem., 45(12), 2615-2623.
        - Daina, A., & Zoete, V. (2016). ChemMedChem, 11(11), 1117-1121.
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
    Calculates all relevant molecular properties for a given molecule.

    Scientific Rationale:
        This function aggregates all property calculations into a single call, providing a comprehensive
        profile of the molecule for downstream scientific analysis, scoring, or database storage.

    Args:
        mol_data (str): Molecular data as a string (SMILES, MOL, SDF).
        input_format (str): Format of the input data ('smiles', 'mol', 'sdf').

    Returns:
        Dict[str, Any]: Dictionary of all calculated properties, including hydrogen bonding, logP, TPSA,
                        molecular properties, functional groups, permeability, and structure identifiers.

    Notes:
        - Returns error information if parsing fails.
        - Adds SMILES, InChI, and InChIKey for structure reference.

    References:
        - See individual property functions for literature references.
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
    Generates an SVG image of a molecule for visualization.

    Scientific Rationale:
        Molecular visualization is essential for interpreting chemical structure, functional groups, and
        substructure matches. Highlighting atoms can be used to emphasize reactive sites or substructure matches.

    Args:
        mol_data (str): Molecular data as a string (SMILES, MOL, SDF).
        input_format (str): Format of the input data ('smiles', 'mol', 'sdf').
        width (int): Image width in pixels.
        height (int): Image height in pixels.
        highlight_atoms (List[int], optional): List of atom indices to highlight.

    Returns:
        str: SVG string representation of the molecule.

    Notes:
        - Hydrogens are removed for cleaner visualization.
        - Stereo annotations are included for clarity.
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
    Performs a substructure search to find occurrences of a query pattern within a target molecule.

    Scientific Rationale:
        Substructure searching is fundamental in cheminformatics for identifying motifs, pharmacophores,
        or functional groups within larger molecules. SMARTS patterns allow for flexible substructure queries.

    Args:
        query_mol_data (str): Query molecule data (usually SMARTS pattern).
        target_mol_data (str): Target molecule data to search in.
        query_format (str): Format of the query data ('smarts', 'smiles').
        target_format (str): Format of the target data ('smiles', 'mol', 'sdf').

    Returns:
        Dict[str, Any]: Dictionary with search results, including match status, count, atom indices, and visualization.

    Notes:
        - Returns a visualization with highlighted matches if found.
        - Useful for SAR (structure-activity relationship) studies.

    References:
        - SMARTS language: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
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
    Calculates molecular similarity between two molecules using various fingerprinting methods.

    Scientific Rationale:
        Molecular similarity is a core concept in cheminformatics, used for virtual screening, clustering,
        and SAR analysis. Fingerprints (e.g., Morgan, MACCS) encode molecular structure for rapid comparison.
        Tanimoto and Dice coefficients are standard similarity metrics.

    Args:
        mol1_data (str): First molecule data.
        mol2_data (str): Second molecule data.
        mol1_format (str): Format of the first molecule data.
        mol2_format (str): Format of the second molecule data.
        fingerprint_type (str): Type of fingerprint to use ('morgan', 'maccs', 'topological').

    Returns:
        Dict[str, float]: Dictionary with similarity metrics (Tanimoto, Dice) and fingerprint type.

    References:
        - Rogers, D., & Hahn, M. (2010). Extended-connectivity fingerprints. J. Chem. Inf. Model., 50(5), 742-754.
        - Willett, P. (2014). Chemical similarity searching. J. Chem. Inf. Model., 54(1), 1-2.
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