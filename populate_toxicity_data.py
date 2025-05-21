#!/usr/bin/env python3
"""
Populate toxicity data from Tox21.

This script:
1. Downloads and processes Tox21 assay and chemical data
2. Maps compounds to existing molecules in the database
3. Imports toxicity data and calculates scores
"""

import os
import sys
import logging
import argparse
from supabase import create_client
from config import Config
from chemical_data.toxicity.tox21_client import Tox21Client
from chemical_data.toxicity.toxicity_scorer import ToxicityScorer

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Populate toxicity data from Tox21')
    parser.add_argument('--limit', type=int, default=None, help='Limit the number of compounds to process')
    parser.add_argument('--force-refresh', action='store_true', help='Force refresh of cached data')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run without modifying the database')
    return parser.parse_args()

def main():
    """Main function to populate toxicity data."""
    args = parse_arguments()
    
    # Initialize Supabase client
    config = Config()
    supabase_url = os.environ.get("SUPABASE_URL") or config.SUPABASE_URL
    supabase_key = os.environ.get("SUPABASE_KEY") or config.SUPABASE_KEY
    
    if not supabase_url or not supabase_key:
        logger.error("Supabase URL or key not found")
        return False
    
    supabase = create_client(supabase_url, supabase_key)
    
    # Initialize Tox21 client and toxicity scorer
    tox21_client = Tox21Client(supabase)
    toxicity_scorer = ToxicityScorer(supabase)
    
    # Download and import Tox21 assay data
    logger.info("Downloading and importing Tox21 assay data")
    if not args.dry_run:
        assay_file = tox21_client.download_assay_data(force_refresh=args.force_refresh)
        assays_imported = tox21_client.import_assays(assay_file)
        logger.info(f"Imported {assays_imported} Tox21 assays")
    else:
        logger.info("Dry run - would import Tox21 assay data")
    
    # Download and import Tox21 chemical data
    logger.info("Downloading and importing Tox21 chemical data")
    if not args.dry_run:
        chemical_file = tox21_client.download_chemical_data(force_refresh=args.force_refresh)
        stats = tox21_client.import_chemical_data(chemical_file, max_compounds=args.limit)
        logger.info(f"Chemical data import stats: {stats}")
    else:
        logger.info("Dry run - would import Tox21 chemical data")
    
    # Calculate toxicity scores
    logger.info("Calculating toxicity scores for molecules")
    if not args.dry_run:
        scores_calculated = toxicity_scorer.calculate_toxicity_scores()
        logger.info(f"Calculated toxicity scores for {scores_calculated} molecules")
    else:
        logger.info("Dry run - would calculate toxicity scores")
    
    logger.info("Toxicity data population completed successfully")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)