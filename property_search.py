#!/usr/bin/env python3
"""
Property-based molecule search script for CryoProtect
This script provides a convenient command-line interface for searching molecules
by their properties using RDKit.
"""

import os
import sys
import argparse
import logging
import json
from typing import Dict, List, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Flag to track if we're using real or mock RDKit
RDKIT_AVAILABLE = False

def setup_rdkit():
    """Set up RDKit, falling back to mock if needed"""
    global RDKIT_AVAILABLE
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski
        logger.info("Using real RDKit")
        RDKIT_AVAILABLE = True
    except ImportError:
        logger.warning("RDKit not available, trying mock implementation...")
        try:
            # Try to import the enhanced mock RDKit
            current_dir = os.path.dirname(os.path.abspath(__file__))
            mock_path = os.path.join(current_dir, "enhanced_mock_rdkit.py")
            
            if os.path.exists(mock_path):
                sys.path.insert(0, current_dir)
                import enhanced_mock_rdkit
                mock_dir = enhanced_mock_rdkit.create_mock_rdkit()
                
                # Ensure mock dir is at the beginning of Python path
                if mock_dir in sys.path:
                    sys.path.remove(mock_dir)
                sys.path.insert(0, mock_dir)
                
                # Clear any cached imports
                if 'rdkit' in sys.modules:
                    del sys.modules['rdkit']
                for name in list(sys.modules.keys()):
                    if name.startswith('rdkit.'):
                        del sys.modules[name]
                
                # Import the mock RDKit modules
                from rdkit import Chem
                from rdkit.Chem import Descriptors, Lipinski
                
                logger.info("Using enhanced mock RDKit")
                RDKIT_AVAILABLE = False
            else:
                logger.error("Enhanced mock RDKit not found at %s", mock_path)
                RDKIT_AVAILABLE = False
                return False
        except Exception as e:
            logger.error("Failed to set up mock RDKit: %s", str(e))
            RDKIT_AVAILABLE = False
            return False
    
    return True

def load_molecules(input_file: str) -> List[Dict[str, Any]]:
    """Load molecules from a JSON file"""
    try:
        with open(input_file, 'r') as f:
            data = json.load(f)
        
        # Check if the input is a list or a dict
        if isinstance(data, list):
            return data
        elif isinstance(data, dict) and 'molecules' in data:
            return data['molecules']
        else:
            logger.error("Invalid JSON format: expected list of molecules or dict with 'molecules' key")
            return []
    except Exception as e:
        logger.error("Error loading molecules from file: %s", str(e))
        return []

def parse_smiles_molecules(smiles_list: List[str]) -> List[Dict[str, Any]]:
    """Parse a list of SMILES strings into molecule dictionaries"""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    
    molecules = []
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("Failed to parse SMILES: %s", smiles)
            continue
        
        # Calculate basic properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        
        molecules.append({
            "name": f"Molecule_{i+1}",
            "smiles": smiles,
            "mol": mol,
            "molecular_weight": mw,
            "logp": logp,
            "h_donors": h_donors,
            "h_acceptors": h_acceptors,
            "h_bond_count": h_donors + h_acceptors
        })
    
    return molecules

def search_molecules(molecules: List[Dict[str, Any]], property_filters: Dict[str, Dict[str, float]]) -> List[Dict[str, Any]]:
    """
    Search molecules based on property filters
    
    Args:
        molecules: List of molecule dictionaries
        property_filters: Dictionary of property filters, e.g.
            {
                "molecular_weight": {"min": 100, "max": 500},
                "logp": {"max": 5}
            }
    
    Returns:
        List of molecules that match all filters
    """
    if not property_filters:
        return molecules
    
    results = []
    
    for mol in molecules:
        matches_all = True
        
        for prop, filter_range in property_filters.items():
            if prop not in mol:
                matches_all = False
                break
            
            value = mol[prop]
            
            # Check minimum value
            if "min" in filter_range and value < filter_range["min"]:
                matches_all = False
                break
            
            # Check maximum value
            if "max" in filter_range and value > filter_range["max"]:
                matches_all = False
                break
        
        if matches_all:
            results.append(mol)
    
    return results

def search_by_substructure(molecules: List[Dict[str, Any]], query: str, query_type: str = "smarts") -> List[Dict[str, Any]]:
    """
    Search molecules by substructure
    
    Args:
        molecules: List of molecule dictionaries
        query: SMARTS or SMILES query string
        query_type: Type of query ("smarts" or "smiles")
    
    Returns:
        List of molecules that contain the substructure
    """
    from rdkit import Chem
    
    # Parse query molecule
    if query_type.lower() == "smarts":
        query_mol = Chem.MolFromSmarts(query)
    else:
        query_mol = Chem.MolFromSmiles(query)
    
    if query_mol is None:
        logger.error("Failed to parse query: %s", query)
        return []
    
    results = []
    
    for mol_data in molecules:
        if "mol" not in mol_data:
            continue
        
        mol = mol_data["mol"]
        if mol.HasSubstructMatch(query_mol):
            results.append(mol_data)
    
    return results

def display_results(molecules: List[Dict[str, Any]], output_format: str = "text", output_file: Optional[str] = None):
    """
    Display the search results
    
    Args:
        molecules: List of molecule dictionaries
        output_format: Output format ("text", "json", or "csv")
        output_file: Optional output file path
    """
    if not molecules:
        print("No molecules found matching the search criteria")
        return
    
    # Create simplified result list for output
    results = []
    for mol in molecules:
        # Skip RDKit molecule object for output
        result = {k: v for k, v in mol.items() if k != "mol"}
        results.append(result)
    
    # Output based on format
    if output_format == "json":
        if output_file:
            with open(output_file, 'w') as f:
                json.dump({"molecules": results}, f, indent=2)
            print(f"Results saved to {output_file}")
        else:
            print(json.dumps({"molecules": results}, indent=2))
    
    elif output_format == "csv":
        import csv
        
        # Determine all possible fields
        all_fields = set()
        for result in results:
            all_fields.update(result.keys())
        
        # Prioritize important fields
        fields = []
        for field in ["name", "smiles", "molecular_weight", "logp", "h_donors", "h_acceptors", "h_bond_count"]:
            if field in all_fields:
                fields.append(field)
                all_fields.remove(field)
        
        # Add remaining fields
        fields.extend(sorted(all_fields))
        
        if output_file:
            with open(output_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerows(results)
            print(f"Results saved to {output_file}")
        else:
            writer = csv.DictWriter(sys.stdout, fieldnames=fields)
            writer.writeheader()
            writer.writerows(results)
    
    else:  # text format
        print(f"Found {len(results)} molecules matching the search criteria:")
        for i, result in enumerate(results, 1):
            print(f"\nMolecule {i}:")
            for key, value in result.items():
                if key != "smiles":  # Print SMILES last for readability
                    print(f"  {key}: {value}")
            print(f"  smiles: {result.get('smiles', 'N/A')}")

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Property-based molecule search")
    
    # Input options
    input_group = parser.add_argument_group("Input Options")
    input_group.add_argument("--input", "-i", help="Input JSON file with molecules")
    input_group.add_argument("--smiles", "-s", nargs="+", help="SMILES strings to search")
    
    # Property filter options
    prop_group = parser.add_argument_group("Property Filters")
    prop_group.add_argument("--mw-min", type=float, help="Minimum molecular weight")
    prop_group.add_argument("--mw-max", type=float, help="Maximum molecular weight")
    prop_group.add_argument("--logp-min", type=float, help="Minimum LogP value")
    prop_group.add_argument("--logp-max", type=float, help="Maximum LogP value")
    prop_group.add_argument("--hbd-min", type=float, help="Minimum H-bond donors")
    prop_group.add_argument("--hbd-max", type=float, help="Maximum H-bond donors")
    prop_group.add_argument("--hba-min", type=float, help="Minimum H-bond acceptors")
    prop_group.add_argument("--hba-max", type=float, help="Maximum H-bond acceptors")
    prop_group.add_argument("--hb-min", type=float, help="Minimum total H-bonds")
    prop_group.add_argument("--hb-max", type=float, help="Maximum total H-bonds")
    
    # Substructure search options
    substruct_group = parser.add_argument_group("Substructure Search")
    substruct_group.add_argument("--smarts", help="SMARTS pattern for substructure search")
    substruct_group.add_argument("--substructure", help="SMILES pattern for substructure search")
    
    # Output options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument("--output", "-o", help="Output file path")
    output_group.add_argument("--format", "-f", choices=["text", "json", "csv"], default="text",
                             help="Output format (text, json, csv)")
    
    # Example command
    parser.epilog = """
    Example commands:
      # Search by molecular weight and LogP
      property_search.py --input molecules.json --mw-min 100 --mw-max 500 --logp-max 5 --format json --output results.json
      
      # Search by SMARTS pattern
      property_search.py --input molecules.json --smarts "[OH]" --format csv
      
      # Process SMILES strings directly
      property_search.py --smiles "CCO" "CC(=O)O" "c1ccccc1" --hb-min 2
    """
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    
    return parser.parse_args()

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Set up RDKit
    if not setup_rdkit():
        logger.error("Failed to set up RDKit or mock implementation")
        return 1
    
    # Load molecules
    molecules = []
    
    if args.input:
        molecules = load_molecules(args.input)
        if not molecules:
            logger.error("No molecules loaded from input file")
            return 1
    elif args.smiles:
        molecules = parse_smiles_molecules(args.smiles)
        if not molecules:
            logger.error("No valid molecules parsed from SMILES")
            return 1
    else:
        logger.error("No input provided. Please specify either --input or --smiles")
        return 1
    
    logger.info("Loaded %d molecules", len(molecules))
    
    # Build property filters
    property_filters = {}
    
    if args.mw_min is not None or args.mw_max is not None:
        property_filters["molecular_weight"] = {}
        if args.mw_min is not None:
            property_filters["molecular_weight"]["min"] = args.mw_min
        if args.mw_max is not None:
            property_filters["molecular_weight"]["max"] = args.mw_max
    
    if args.logp_min is not None or args.logp_max is not None:
        property_filters["logp"] = {}
        if args.logp_min is not None:
            property_filters["logp"]["min"] = args.logp_min
        if args.logp_max is not None:
            property_filters["logp"]["max"] = args.logp_max
    
    if args.hbd_min is not None or args.hbd_max is not None:
        property_filters["h_donors"] = {}
        if args.hbd_min is not None:
            property_filters["h_donors"]["min"] = args.hbd_min
        if args.hbd_max is not None:
            property_filters["h_donors"]["max"] = args.hbd_max
    
    if args.hba_min is not None or args.hba_max is not None:
        property_filters["h_acceptors"] = {}
        if args.hba_min is not None:
            property_filters["h_acceptors"]["min"] = args.hba_min
        if args.hba_max is not None:
            property_filters["h_acceptors"]["max"] = args.hba_max
    
    if args.hb_min is not None or args.hb_max is not None:
        property_filters["h_bond_count"] = {}
        if args.hb_min is not None:
            property_filters["h_bond_count"]["min"] = args.hb_min
        if args.hb_max is not None:
            property_filters["h_bond_count"]["max"] = args.hb_max
    
    # Apply property filters
    if property_filters:
        logger.info("Applying property filters: %s", property_filters)
        molecules = search_molecules(molecules, property_filters)
        logger.info("%d molecules remain after property filtering", len(molecules))
    
    # Apply substructure search
    if args.smarts:
        logger.info("Applying SMARTS filter: %s", args.smarts)
        molecules = search_by_substructure(molecules, args.smarts, "smarts")
        logger.info("%d molecules remain after SMARTS filtering", len(molecules))
    
    if args.substructure:
        logger.info("Applying substructure filter: %s", args.substructure)
        molecules = search_by_substructure(molecules, args.substructure, "smiles")
        logger.info("%d molecules remain after substructure filtering", len(molecules))
    
    # Display results
    display_results(molecules, args.format, args.output)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())