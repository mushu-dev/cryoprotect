#!/usr/bin/env python3
"""
PubChem import script with direct database connection for CryoProtect.
This script sets up the connection parameters properly and runs the import.
"""

import os
import sys
import argparse
import logging
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'logs/pubchem_import_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description="Import cryoprotectant data from PubChem with direct connection")
    parser.add_argument("--host", required=True, help="Database host")
    parser.add_argument("--port", default="5432", help="Database port")
    parser.add_argument("--user", required=True, help="Database user")
    parser.add_argument("--password", required=True, help="Database password")
    parser.add_argument("--dbname", default="postgres", help="Database name")
    parser.add_argument("--target", type=int, default=100, help="Target number of compounds to import")
    parser.add_argument("--api-delay", type=float, default=1.0, help="Base delay between API calls (seconds)")
    parser.add_argument("--checkpoint", type=str, help="Path to checkpoint file to resume from")
    parser.add_argument("--batch-size", type=int, default=10, help="Number of compounds to process in each batch")
    args = parser.parse_args()
    
    # Set connection environment variables
    os.environ["SUPABASE_DB_HOST"] = args.host
    os.environ["SUPABASE_DB_PORT"] = args.port
    os.environ["SUPABASE_DB_USER"] = args.user
    os.environ["SUPABASE_DB_PASSWORD"] = args.password
    os.environ["SUPABASE_DB_NAME"] = args.dbname
    
    logger.info(f"Set up database connection to {args.host}:{args.port} as user {args.user}")
    
    # Import the actual script
    import simple_pubchem_import
    
    # Set up arguments for the import
    sys.argv = [
        "simple_pubchem_import.py",
        "--target", str(args.target),
        "--api-delay", str(args.api_delay),
        "--batch-size", str(args.batch_size)
    ]
    
    if args.checkpoint:
        sys.argv.extend(["--checkpoint", args.checkpoint])
    
    # Run the import
    simple_pubchem_import.main()

if __name__ == "__main__":
    main()