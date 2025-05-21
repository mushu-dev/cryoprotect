# cryoprotectant_identifiers.py
import json
import os
import sys
import logging
from typing import Dict, List, Any, Optional, Set, Tuple

logger = logging.getLogger(__name__)

MASTER_LIST_PATH = os.path.join('data', 'cryoprotectant_master_list.json')

class CryoprotectantIdentifierManager:
    """
    Manages cryoprotectant identifiers across different data sources.
    Provides mapping between different identifier systems and validation.
    """
    _instance = None
    
    def __init__(self, master_list_path: str = MASTER_LIST_PATH):
        """
        Initialize with path to master identifier list.
        
        Args:
            master_list_path: Path to the JSON file containing master identifier list
        """
        self.master_list_path = master_list_path
        self.identifiers = {}
        self.names = {}
        self.pubchem_cids = set()
        self.chembl_ids = set()
        self.cas_numbers = set()
        self.inchi_keys = set()
        self.smiles = set()
        self._load_identifiers()
    
    @classmethod
    def get_instance(cls, master_list_path: str = MASTER_LIST_PATH) -> 'CryoprotectantIdentifierManager':
        """
        Get singleton instance.
        
        Args:
            master_list_path: Path to master identifier list
            
        Returns:
            CryoprotectantIdentifierManager instance
        """
        if cls._instance is None:
            cls._instance = cls(master_list_path)
        return cls._instance
    
    def _load_identifiers(self) -> None:
        """Load identifiers from master list file."""
        try:
            if not os.path.exists(self.master_list_path):
                logger.warning("Master list file not found at %s. Creating empty list.", 
                              self.master_list_path)
                self.identifiers = {}
                return
                
            with open(self.master_list_path, 'r') as f:
                self.identifiers = json.load(f)
            
            # Build lookup indexes
            for internal_id, data in self.identifiers.items():
                if 'pubchem_cid' in data and data['pubchem_cid']:
                    self.pubchem_cids.add(str(data['pubchem_cid']))
                
                if 'chembl_id' in data and data['chembl_id']:
                    self.chembl_ids.add(data['chembl_id'])
                
                if 'cas_number' in data and data['cas_number']:
                    self.cas_numbers.add(data['cas_number'])
                
                if 'inchi_key' in data and data['inchi_key']:
                    self.inchi_keys.add(data['inchi_key'])
                
                if 'smiles' in data and data['smiles']:
                    self.smiles.add(data['smiles'])
                
                if 'names' in data and data['names']:
                    for name in data['names']:
                        self.names[name.lower()] = internal_id
            
            logger.info("Loaded %d cryoprotectant identifiers", len(self.identifiers))
            logger.info("Indexed %d PubChem CIDs, %d ChEMBL IDs, %d CAS numbers, %d InChI Keys, %d SMILES",
                       len(self.pubchem_cids), len(self.chembl_ids), len(self.cas_numbers),
                       len(self.inchi_keys), len(self.smiles))
                       
        except Exception as e:
            logger.error("Failed to load identifiers: %s", str(e))
            # Initialize empty
            self.identifiers = {}
    
    def save_identifiers(self) -> None:
        """Save identifiers to master list file."""
        try:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.master_list_path), exist_ok=True)
            
            with open(self.master_list_path, 'w') as f:
                json.dump(self.identifiers, f, indent=2)
            
            logger.info("Saved %d cryoprotectant identifiers to %s", 
                       len(self.identifiers), self.master_list_path)
        except Exception as e:
            logger.error("Failed to save identifiers: %s", str(e))
    
    def add_molecule(self, internal_id: str, data: Dict[str, Any]) -> None:
        """
        Add a molecule to the identifier list.
        
        Args:
            internal_id: Internal identifier for the molecule
            data: Dictionary containing molecule identifier data
        """
        self.identifiers[internal_id] = data
        
        # Update lookup indexes
        if 'pubchem_cid' in data and data['pubchem_cid']:
            self.pubchem_cids.add(str(data['pubchem_cid']))
        
        if 'chembl_id' in data and data['chembl_id']:
            self.chembl_ids.add(data['chembl_id'])
        
        if 'cas_number' in data and data['cas_number']:
            self.cas_numbers.add(data['cas_number'])
        
        if 'inchi_key' in data and data['inchi_key']:
            self.inchi_keys.add(data['inchi_key'])
        
        if 'smiles' in data and data['smiles']:
            self.smiles.add(data['smiles'])
        
        if 'names' in data and data['names']:
            for name in data['names']:
                self.names[name.lower()] = internal_id
    
    def get_molecule_by_internal_id(self, internal_id: str) -> Optional[Dict[str, Any]]:
        """
        Get molecule data by internal identifier.
        
        Args:
            internal_id: Internal identifier for the molecule
            
        Returns:
            Dictionary containing molecule data or None if not found
        """
        return self.identifiers.get(internal_id)
    
    def get_internal_id_by_pubchem_cid(self, pubchem_cid: str) -> Optional[str]:
        """
        Get internal identifier by PubChem CID.
        
        Args:
            pubchem_cid: PubChem CID
            
        Returns:
            Internal identifier or None if not found
        """
        for internal_id, data in self.identifiers.items():
            if 'pubchem_cid' in data and str(data['pubchem_cid']) == str(pubchem_cid):
                return internal_id
        return None
    
    def get_internal_id_by_chembl_id(self, chembl_id: str) -> Optional[str]:
        """
        Get internal identifier by ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            Internal identifier or None if not found
        """
        for internal_id, data in self.identifiers.items():
            if 'chembl_id' in data and data['chembl_id'] == chembl_id:
                return internal_id
        return None
    
    def get_internal_id_by_name(self, name: str) -> Optional[str]:
        """
        Get internal identifier by molecule name.
        
        Args:
            name: Molecule name
            
        Returns:
            Internal identifier or None if not found
        """
        return self.names.get(name.lower())
    
    def get_internal_id_by_inchi_key(self, inchi_key: str) -> Optional[str]:
        """
        Get internal identifier by InChI Key.
        
        Args:
            inchi_key: InChI Key
            
        Returns:
            Internal identifier or None if not found
        """
        for internal_id, data in self.identifiers.items():
            if 'inchi_key' in data and data['inchi_key'] == inchi_key:
                return internal_id
        return None
    
    def resolve_identifier(self, 
                          pubchem_cid: Optional[str] = None,
                          chembl_id: Optional[str] = None,
                          name: Optional[str] = None,
                          cas_number: Optional[str] = None, 
                          inchi_key: Optional[str] = None,
                          smiles: Optional[str] = None) -> Tuple[Optional[str], float]:
        """
        Resolve molecule identifier from any available identifier.
        
        Args:
            pubchem_cid: PubChem CID
            chembl_id: ChEMBL ID
            name: Molecule name
            cas_number: CAS Registry Number
            inchi_key: InChI Key
            smiles: SMILES string
            
        Returns:
            Tuple of (internal_id, confidence_score) or (None, 0.0) if not found
        """
        # Check exact matches
        if pubchem_cid:
            internal_id = self.get_internal_id_by_pubchem_cid(pubchem_cid)
            if internal_id:
                return internal_id, 1.0
        
        if chembl_id:
            internal_id = self.get_internal_id_by_chembl_id(chembl_id)
            if internal_id:
                return internal_id, 1.0
        
        if inchi_key:
            internal_id = self.get_internal_id_by_inchi_key(inchi_key)
            if internal_id:
                return internal_id, 1.0
        
        if name:
            internal_id = self.get_internal_id_by_name(name)
            if internal_id:
                return internal_id, 0.9  # Names can be ambiguous
        
        # No match found
        return None, 0.0
    
    def get_all_pubchem_cids(self) -> List[str]:
        """
        Get list of all PubChem CIDs.
        
        Returns:
            List of PubChem CIDs
        """
        return list(self.pubchem_cids)
    
    def get_all_chembl_ids(self) -> List[str]:
        """
        Get list of all ChEMBL IDs.
        
        Returns:
            List of ChEMBL IDs
        """
        return list(self.chembl_ids)

def initialize_cryoprotectant_list() -> None:
    """Initialize the cryoprotectant list with commonly used molecules."""
    manager = CryoprotectantIdentifierManager.get_instance()
    
    # Core cryoprotectants
    common_cryoprotectants = [
        {
            "internal_id": "CRYO001",
            "pubchem_cid": "962",
            "chembl_id": "CHEMBL388978",
            "cas_number": "56-81-5",
            "names": ["Glycerol", "Glycerin", "1,2,3-Propanetriol"],
            "inchi_key": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
            "smiles": "C(C(CO)O)O",
            "formula": "C3H8O3",
            "molecular_weight": 92.09,
            "category": "polyol"
        },
        {
            "internal_id": "CRYO002", 
            "pubchem_cid": "5988", 
            "chembl_id": "CHEMBL1098659", 
            "cas_number": "67-68-5", 
            "names": ["Dimethyl sulfoxide", "DMSO", "Methyl sulfoxide"], 
            "inchi_key": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N", 
            "smiles": "CS(=O)C", 
            "formula": "C2H6OS", 
            "molecular_weight": 78.13, 
            "category": "organosulfur"
        },
        {
            "internal_id": "CRYO003",
            "pubchem_cid": "6342",
            "chembl_id": "CHEMBL66195",
            "cas_number": "107-95-9",
            "names": ["beta-Alanine", "3-Aminopropanoic acid", "3-Aminopropionic acid"],
            "inchi_key": "UCMIRNVEIXFBKS-UHFFFAOYSA-N",
            "smiles": "C(CN)C(=O)O",
            "formula": "C3H7NO2",
            "molecular_weight": 89.09,
            "category": "amino acid"
        },
        {
            "internal_id": "CRYO004",
            "pubchem_cid": "1030",
            "chembl_id": "CHEMBL500033",
            "cas_number": "75-65-0",
            "names": ["tert-Butanol", "t-Butyl alcohol", "2-Methyl-2-propanol"],
            "inchi_key": "DKGAVHZHDRPRBM-UHFFFAOYSA-N",
            "smiles": "CC(C)(C)O",
            "formula": "C4H10O",
            "molecular_weight": 74.12,
            "category": "alcohol"
        },
        {
            "internal_id": "CRYO005",
            "pubchem_cid": "6057",
            "chembl_id": "CHEMBL1487",
            "cas_number": "57-13-6",
            "names": ["Urea", "Carbamide", "Carbonyl diamide"],
            "inchi_key": "XSQUKJJJFZCRTK-UHFFFAOYSA-N",
            "smiles": "C(=O)(N)N",
            "formula": "CH4N2O",
            "molecular_weight": 60.06,
            "category": "amide"
        }
    ]
    
    # Add molecules to manager
    for molecule in common_cryoprotectants:
        manager.add_molecule(molecule["internal_id"], molecule)
    
    # Save to file
    manager.save_identifiers()
    
    logger.info("Initialized cryoprotectant list with %d molecules", 
               len(common_cryoprotectants))

def validate_identifiers():
    """Validate identifier consistency across different systems."""
    logger.info("Validating cryoprotectant identifiers...")
    
    # Get the identifier manager
    manager = CryoprotectantIdentifierManager.get_instance()
    
    # Check if the manager has loaded identifiers
    if not manager.identifiers:
        logger.warning("No identifiers loaded in CryoprotectantIdentifierManager")
        return False
    
    # Validation results
    results = {
        "total_molecules": len(manager.identifiers),
        "pubchem_ids": len(manager.pubchem_cids),
        "chembl_ids": len(manager.chembl_ids),
        "cas_numbers": len(manager.cas_numbers),
        "inchi_keys": len(manager.inchi_keys),
        "smiles": len(manager.smiles),
        "names": len(manager.names),
        "cross_reference_issues": []
    }
    
    # Check for cross-reference consistency
    for internal_id, data in manager.identifiers.items():
        issues = []
        
        # Check PubChem CID
        if 'pubchem_cid' in data and data['pubchem_cid']:
            pubchem_cid = str(data['pubchem_cid'])
            resolved_id, confidence = manager.resolve_identifier(pubchem_cid=pubchem_cid)
            if resolved_id != internal_id:
                issues.append(f"PubChem CID {pubchem_cid} resolves to {resolved_id} instead of {internal_id}")
        
        # Check ChEMBL ID
        if 'chembl_id' in data and data['chembl_id']:
            chembl_id = data['chembl_id']
            resolved_id, confidence = manager.resolve_identifier(chembl_id=chembl_id)
            if resolved_id != internal_id:
                issues.append(f"ChEMBL ID {chembl_id} resolves to {resolved_id} instead of {internal_id}")
        
        # Check InChI Key
        if 'inchi_key' in data and data['inchi_key']:
            inchi_key = data['inchi_key']
            resolved_id, confidence = manager.resolve_identifier(inchi_key=inchi_key)
            if resolved_id != internal_id:
                issues.append(f"InChI Key {inchi_key} resolves to {resolved_id} instead of {internal_id}")
        
        # If there are issues, add to results
        if issues:
            results["cross_reference_issues"].append({
                "internal_id": internal_id,
                "name": data.get('names', ['Unknown'])[0],
                "issues": issues
            })
    
    # Check for database consistency if possible
    try:
        import psycopg2
        from psycopg2.extras import RealDictCursor
        import os
        
        # Try to connect to database
        db_params = {
            'host': os.environ.get('SUPABASE_DB_HOST'),
            'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
            'user': os.environ.get('SUPABASE_DB_USER'),
            'password': os.environ.get('SUPABASE_DB_PASSWORD'),
            'dbname': os.environ.get('SUPABASE_DB_NAME')
        }
        
        # Check if all required parameters are present
        if all(db_params.values()):
            conn = psycopg2.connect(**db_params)
            cursor = conn.cursor(cursor_factory=RealDictCursor)
            
            # Check database schema to determine how identifiers are stored
            cursor.execute("SELECT column_name FROM information_schema.columns WHERE table_name = 'molecules'")
            columns = [row['column_name'] for row in cursor.fetchall()]
            
            # Determine if identifiers are stored in columns or as JSON properties
            has_pubchem_column = "pubchem_cid" in columns
            has_chembl_column = "chembl_id" in columns
            has_properties_column = "properties" in columns
            
            # Count molecules with both PubChem and ChEMBL identifiers
            if has_pubchem_column and has_chembl_column:
                # Identifiers are stored in separate columns
                cursor.execute("""
                SELECT
                  count(*) as total_molecules,
                  sum(CASE WHEN pubchem_cid IS NOT NULL AND chembl_id IS NOT NULL THEN 1 ELSE 0 END) as molecules_with_both_sources,
                  sum(CASE WHEN pubchem_cid IS NOT NULL AND chembl_id IS NULL THEN 1 ELSE 0 END) as pubchem_only,
                  sum(CASE WHEN pubchem_cid IS NULL AND chembl_id IS NOT NULL THEN 1 ELSE 0 END) as chembl_only
                FROM molecules
                """)
            elif has_properties_column:
                # Identifiers are stored in JSON properties column
                cursor.execute("""
                SELECT
                  count(*) as total_molecules,
                  sum(CASE WHEN properties->>'pubchem' IS NOT NULL AND properties->>'chembl' IS NOT NULL THEN 1 ELSE 0 END) as molecules_with_both_sources,
                  sum(CASE WHEN properties->>'pubchem' IS NOT NULL AND properties->>'chembl' IS NULL THEN 1 ELSE 0 END) as pubchem_only,
                  sum(CASE WHEN properties->>'pubchem' IS NULL AND properties->>'chembl' IS NOT NULL THEN 1 ELSE 0 END) as chembl_only
                FROM molecules
                """)
            
            row = cursor.fetchone()
            if row:
                results["database"] = {
                    "total_molecules": row["total_molecules"],
                    "molecules_with_both_sources": row["molecules_with_both_sources"],
                    "pubchem_only": row["pubchem_only"],
                    "chembl_only": row["chembl_only"]
                }
            
            cursor.close()
            conn.close()
    except Exception as e:
        logger.warning(f"Could not check database consistency: {e}")
    
    # Print validation results
    print("\n" + "=" * 60)
    print("Cryoprotectant Identifier Validation Results")
    print("=" * 60)
    print(f"Total molecules: {results['total_molecules']}")
    print(f"PubChem CIDs: {results['pubchem_ids']}")
    print(f"ChEMBL IDs: {results['chembl_ids']}")
    print(f"CAS Numbers: {results['cas_numbers']}")
    print(f"InChI Keys: {results['inchi_keys']}")
    print(f"SMILES: {results['smiles']}")
    print(f"Names: {results['names']}")
    
    if "database" in results:
        print("\nDatabase Statistics:")
        print(f"Total molecules in database: {results['database']['total_molecules']}")
        print(f"Molecules with both PubChem and ChEMBL IDs: {results['database']['molecules_with_both_sources']}")
        print(f"Molecules with PubChem ID only: {results['database']['pubchem_only']}")
        print(f"Molecules with ChEMBL ID only: {results['database']['chembl_only']}")
    
    if results["cross_reference_issues"]:
        print("\nCross-Reference Issues:")
        for issue in results["cross_reference_issues"]:
            print(f"  {issue['name']} ({issue['internal_id']}):")
            for problem in issue['issues']:
                print(f"    - {problem}")
    else:
        print("\nNo cross-reference issues found.")
    
    print("=" * 60)
    
    # Save results to file
    os.makedirs('reports', exist_ok=True)
    with open('reports/identifier_validation.json', 'w') as f:
        import json
        json.dump(results, f, indent=2)
    
    logger.info(f"Validation results saved to reports/identifier_validation.json")
    
    return len(results["cross_reference_issues"]) == 0

# Command-line interface
if __name__ == "__main__":
    import argparse
    
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Cryoprotectant Identifier Manager")
    parser.add_argument("--validate", action="store_true", help="Validate identifier consistency")
    parser.add_argument("--initialize", action="store_true", help="Initialize the cryoprotectant list")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()
    
    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level)
    logger = logging.getLogger(__name__)
    
    if args.validate:
        # Run validation
        logger.info("Running identifier validation...")
        success = validate_identifiers()
        sys.exit(0 if success else 1)
    else:
        # Default: initialize the list
        logger.info("Initializing cryoprotectant identifier list...")
        initialize_cryoprotectant_list()
        logger.info("Done!")