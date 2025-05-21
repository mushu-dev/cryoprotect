# test_cryoprotectant_identifiers.py
import logging
from cryoprotectant_identifiers import CryoprotectantIdentifierManager, initialize_cryoprotectant_list

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_identifier_manager():
    """Test the cryoprotectant identifier manager."""
    # Initialize with common cryoprotectants
    initialize_cryoprotectant_list()
    
    # Get instance
    manager = CryoprotectantIdentifierManager.get_instance()
    
    # Test lookup by different identifiers
    glycerol_by_cid = manager.get_internal_id_by_pubchem_cid("962")
    glycerol_by_name = manager.get_internal_id_by_name("Glycerol")
    glycerol_by_inchi = manager.get_internal_id_by_inchi_key("PEDCQBHIVMGVHV-UHFFFAOYSA-N")
    
    logger.info(f"Glycerol by CID: {glycerol_by_cid}")
    logger.info(f"Glycerol by name: {glycerol_by_name}")
    logger.info(f"Glycerol by InChI key: {glycerol_by_inchi}")
    
    # Test resolve_identifier
    dmso_id, confidence = manager.resolve_identifier(pubchem_cid="5988")
    logger.info(f"DMSO resolved to {dmso_id} with confidence {confidence}")
    
    # Test adding a new molecule
    manager.add_molecule("CRYO006", {
        "internal_id": "CRYO006",
        "pubchem_cid": "702",
        "chembl_id": "CHEMBL388978",
        "cas_number": "50-70-4",
        "names": ["Sorbitol", "D-Glucitol", "D-Sorbitol"],
        "inchi_key": "FBPFZTCFMRRESA-JGWLITMVSA-N",
        "smiles": "C(C(C(C(C(CO)O)O)O)O)O",
        "formula": "C6H14O6",
        "molecular_weight": 182.17,
        "category": "polyol"
    })
    
    # Verify the new molecule
    sorbitol_id = manager.get_internal_id_by_name("Sorbitol")
    logger.info(f"Sorbitol ID: {sorbitol_id}")
    
    # Save the updated list
    manager.save_identifiers()
    
    return True

if __name__ == "__main__":
    logger.info("Testing cryoprotectant identifier manager...")
    test_identifier_manager()
    logger.info("Cryoprotectant identifier test completed!")