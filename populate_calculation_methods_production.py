"""
populate_calculation_methods_production.py

Populates the 'calculation_method' table with a curated, production-grade list of experimental, computational, and machine learning methods for CryoProtect v2.

Approach:
- Connects to the Supabase database using the client.
- Defines a curated list of calculation methods, each with name, description, method_type, and reference to literature or documentation.
- Inserts the methods into the calculation_method table (or outputs data if dry_run).
- Documents sources and process for reproducibility.

Requirements:
- supabase-py (install via pip if needed)
- Python 3.7+

Usage:
    python populate_calculation_methods_production.py [--dry-run]

Authoritative sources:
- Fahy, G. M., et al. "Vitrification as an approach to cryopreservation." Cryobiology 21.4 (1984): 407-426.
- Rall, W. F., & Fahy, G. M. "Ice-free cryopreservation of mouse embryos at −196°C by vitrification." Nature 313.6003 (1985): 573-575.
- PubChem, ChemSpider, RDKit documentation.
- Scikit-learn documentation for ML methods.
- Additional protocols from Nature Protocols, Cryobiology, and peer-reviewed literature.

This script is suitable for direct production use and is fully reproducible.
"""

import argparse
import uuid
import os
import logging
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client, Client
from service_role_helper import get_supabase_client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Curated list of calculation methods
CALCULATION_METHODS = [
    {
        "name": "Differential Scanning Calorimetry (DSC)",
        "description": "Experimental method for measuring glass transition and melting points of cryoprotectant solutions.",
        "method_type": "experimental",
        "reference": "Fahy et al., Cryobiology 21.4 (1984): 407-426; Kasai et al., Biology of Reproduction 28.3 (1983): 687-693."
    },
    {
        "name": "Manual Freezing Point Depression",
        "description": "Experimental measurement of freezing point by controlled cooling and observation.",
        "method_type": "experimental",
        "reference": "Rall & Fahy, Nature 313 (1985): 573-575."
    },
    {
        "name": "PubChem Data Extraction",
        "description": "Computational extraction of molecular properties from PubChem database.",
        "method_type": "computational",
        "reference": "PubChem (https://pubchem.ncbi.nlm.nih.gov)"
    },
    {
        "name": "RDKit Property Calculation",
        "description": "Computational calculation of molecular descriptors (LogP, TPSA, etc.) using RDKit.",
        "method_type": "computational",
        "reference": "RDKit documentation (https://www.rdkit.org/docs/)"
    },
    {
        "name": "Random Forest Regression (ML)",
        "description": "Machine learning model for predicting cryoprotectant properties using molecular descriptors.",
        "method_type": "ML",
        "reference": "scikit-learn documentation (https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html)"
    },
    {
        "name": "Gradient Boosting Regression (ML)",
        "description": "Machine learning model for predicting cryoprotectant properties using molecular descriptors.",
        "method_type": "ML",
        "reference": "scikit-learn documentation (https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.GradientBoostingRegressor.html)"
    },
    {
        "name": "Neural Network Regression (ML)",
        "description": "Machine learning model for predicting cryoprotectant properties using molecular descriptors.",
        "method_type": "ML",
        "reference": "scikit-learn documentation (https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPRegressor.html)"
    }
]

def insert_calculation_method(supabase, method, created_by=None, dry_run=False):
    """Insert a calculation method into the database using Supabase client."""
    method_id = str(uuid.uuid4())
    now = datetime.utcnow()
    
    # Prepare method data
    method_data = {
        "id": method_id,
        "name": method["name"],
        "description": method["description"],
        "method_type": method["method_type"],
        "reference": method["reference"],
        "created_by": created_by,
        "created_at": now.isoformat(),
        "updated_at": now.isoformat()
    }
    
    if not dry_run:
        # Insert method
        response = supabase.table("calculation_method").insert(method_data).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error inserting calculation method {method['name']}: {response.error}")
            return False
        logger.info(f"Inserted calculation method {method['name']} with ID: {method_id}")
        return True
    else:
        logger.info(f"DRY RUN: Would insert calculation method {method['name']} with ID: {method_id}")
        return True

def main():
    parser = argparse.ArgumentParser(description="Populate calculation_method table with production-grade methods.")
    parser.add_argument("--dry-run", action="store_true", help="Print data instead of executing")
    parser.add_argument("--created-by", help="UUID of the user creating these records (optional)")
    args = parser.parse_args()

    # Connect to Supabase
    try:
        supabase = get_supabase_client()
        logger.info("Connected to Supabase using service role")
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        # Fall back to direct connection
        supabase_url = os.getenv("SUPABASE_URL")
        supabase_key = os.getenv("SUPABASE_KEY")
        if not supabase_url or not supabase_key:
            logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
            return
        supabase = create_client(supabase_url, supabase_key)
        logger.info("Connected to Supabase directly")

    try:
        success_count = 0
        failure_count = 0
        
        for method in CALCULATION_METHODS:
            if insert_calculation_method(supabase, method, created_by=args.created_by, dry_run=args.dry_run):
                success_count += 1
            else:
                failure_count += 1
        
        if not args.dry_run:
            logger.info(f"Calculation method population complete. Successes: {success_count}, Failures: {failure_count}")
        else:
            logger.info(f"Dry run complete. No changes made to the database. Would have processed: {success_count} calculation methods")
    except Exception as e:
        logger.error(f"Error in calculation method population: {str(e)}")

if __name__ == "__main__":
    main()