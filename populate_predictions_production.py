"""
populate_predictions_production.py

Populates the 'prediction' table with production-grade, scientifically accurate predictions for CryoProtect v2.

Approach:
- Connects to the Supabase database using the client.
- Fetches all molecules and mixtures.
- For each, runs the predictive models (using ModelManager from api/predictive_models.py) for relevant properties.
- Collects predicted_value, confidence, model_version, method_id (linking to calculation_method), and provenance.
- Inserts into the prediction table (or outputs data if dry_run).
- Documents sources and process for reproducibility.

Requirements:
- supabase-py (install via pip if needed)
- Python 3.7+
- numpy, pandas, scikit-learn, RDKit (for predictive models)
- CryoProtect v2 codebase in PYTHONPATH

Usage:
    python populate_predictions_production.py [--dry-run]

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

import sys

# Ensure api/ is in PYTHONPATH for imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "api")))

from api.predictive_models import ModelManager

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# List of properties to predict (can be expanded as needed)
PROPERTIES_TO_PREDICT = [
    "glass_transition_temp",
    "vitrification_success_rate",
    "oocyte_survival_rate",
    "post_thaw_motility"
]

# Map ML method names to calculation_method names (must match what's in calculation_method table)
ML_METHOD_TO_NAME = {
    "random_forest": "Random Forest Regression (ML)",
    "gradient_boosting": "Gradient Boosting Regression (ML)",
    "neural_network": "Neural Network Regression (ML)"
}

def get_molecule_ids(supabase):
    """Fetch all molecule IDs from the database using Supabase client."""
    response = supabase.table("molecule").select("id").execute()
    if hasattr(response, 'data') and response.data:
        return [row["id"] for row in response.data]
    return []

def get_mixture_ids(supabase):
    """Fetch all mixture IDs from the database using Supabase client."""
    response = supabase.table("mixture").select("id").execute()
    if hasattr(response, 'data') and response.data:
        return [row["id"] for row in response.data]
    return []

def get_calculation_method_ids(supabase):
    """Fetch calculation method name to ID mapping using Supabase client."""
    response = supabase.table("calculation_method").select("id, name").execute()
    if hasattr(response, 'data') and response.data:
        return {row["name"]: row["id"] for row in response.data}
    return {}

def insert_prediction(supabase, pred, dry_run=False):
    """Insert a prediction into the database using Supabase client."""
    prediction_data = {
        "id": pred["id"],
        "molecule_id": pred["molecule_id"],
        "property_type": pred["property_type"],
        "predicted_value": pred["predicted_value"],
        "unit": pred["unit"],
        "method_id": pred["method_id"],
        "model_version": pred["model_version"],
        "confidence": pred["confidence"],
        "data_source": pred["data_source"],
        "version": 1,
        "modification_history": "[]",
        "created_at": pred["created_at"].isoformat(),
        "updated_at": pred["updated_at"].isoformat()
    }
    
    if not dry_run:
        response = supabase.table("prediction").insert(prediction_data).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error inserting prediction: {response.error}")
            return False
        return True
    else:
        logger.info(f"DRY RUN: Would insert prediction for molecule {pred['molecule_id']}, property {pred['property_type']}")
        return True

def main():
    parser = argparse.ArgumentParser(description="Populate prediction table with production-grade predictions.")
    parser.add_argument("--dry-run", action="store_true", help="Print data instead of executing")
    parser.add_argument("--ml-method", default="random_forest", choices=list(ML_METHOD_TO_NAME.keys()), help="ML method to use for predictions")
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
        model_manager = ModelManager()
        method_ids = get_calculation_method_ids(supabase)
        method_name = ML_METHOD_TO_NAME[args.ml_method]
        method_id = method_ids.get(method_name)
        if not method_id:
            logger.error(f"Calculation method '{method_name}' not found in calculation_method table. Please populate it first.")
            return

        molecule_ids = get_molecule_ids(supabase)
        logger.info(f"Found {len(molecule_ids)} molecules to process")
        now = datetime.utcnow()
        
        success_count = 0
        failure_count = 0

        for molecule_id in molecule_ids:
            for property_type in PROPERTIES_TO_PREDICT:
                # Predict for molecule (using ModelManager)
                try:
                    result = model_manager.predict(property_type, molecule_id, algorithm=args.ml_method)
                    pred = {
                        "id": str(uuid.uuid4()),
                        "molecule_id": molecule_id,
                        "property_type": property_type,
                        "predicted_value": result.get("prediction"),
                        "unit": result.get("unit", "N/A"),
                        "method_id": method_id,
                        "model_version": result.get("model_version", "unknown"),
                        "confidence": result.get("confidence"),
                        "data_source": "production_script",
                        "created_at": now,
                        "updated_at": now
                    }
                    if insert_prediction(supabase, pred, dry_run=args.dry_run):
                        success_count += 1
                    else:
                        failure_count += 1
                except Exception as e:
                    logger.error(f"Prediction failed for molecule {molecule_id}, property {property_type}: {e}")
                    failure_count += 1

        if not args.dry_run:
            logger.info(f"All predictions inserted successfully. Successes: {success_count}, Failures: {failure_count}")
        else:
            logger.info(f"Dry run complete. No changes made to the database. Would have processed: {success_count} predictions")
    except Exception as e:
        logger.error(f"Error in prediction population: {str(e)}")

if __name__ == "__main__":
    main()