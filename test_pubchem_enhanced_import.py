# test_pubchem_enhanced_import.py
import os
import logging
import argparse
from PubChem_CryoProtectants_Supabase_Enhanced import PubChemImporter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description='Test PubChem enhanced import')
    parser.add_argument('--cids', type=str, default="962,5988,6342,1030,6057",
                       help='Comma-separated list of CIDs to test import')
    parser.add_argument('--checkpoint', type=str, 
                       default="checkpoints/pubchem_import_test.json",
                       help='Checkpoint file path')
    parser.add_argument('--report', type=str,
                       default="reports/pubchem_import_test_report.json",
                       help='Report file path')
    args = parser.parse_args()
    
    # Parse CIDs
    cids = [cid.strip() for cid in args.cids.split(',') if cid.strip()]
    
    if not cids:
        logger.error("No CIDs provided")
        return 1
    
    logger.info("Testing import with %d CIDs: %s", len(cids), ', '.join(cids))
    
    # Create importer
    importer = PubChemImporter(
        checkpoint_file=args.checkpoint,
        batch_size=2  # Small batch size for testing
    )
    
    # Import compounds
    importer.import_compounds(cids)
    
    # Generate report
    importer.generate_report(args.report)
    
    return 0

if __name__ == "__main__":
    main()