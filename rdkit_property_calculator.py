#!/usr/bin/env python3
"""
RDKit Property Calculator for CryoProtect

This script calculates molecular properties for cryoprotectants using RDKit and
stores them in the Supabase database. It uses the rdkit_wrapper to ensure
consistent property calculation across environments.

Usage:
  python rdkit_property_calculator.py [options]

Options:
  --known     Process only known cryoprotectants
  --sample N  Process a random sample of N molecules
  --limit N   Limit processing to N molecules
  --help      Show this help message
"""

import os
import sys
import argparse
import uuid
import json
import logging
from datetime import datetime
from dotenv import load_dotenv
from tqdm import tqdm
import supabase
from supabase import create_client, Client

# Import rdkit_wrapper (our unified RDKit interface)
from rdkit_wrapper import (
    create_molecule_from_smiles,
    calculate_properties,
    get_rdkit_status,
    RDKIT_AVAILABLE
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("rdkit_property_calculation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# List of cryoprotectant-specific properties to calculate
CRYOPROTECTANT_PROPERTIES = [
    {
        "name": "h_bond_donor_acceptor_ratio",
        "description": "Ratio of H-bond donors to acceptors",
        "calculation": lambda props: (
            props.get("h_donors", 0) / props.get("h_acceptors", 1)
            if props.get("h_acceptors", 0) > 0 else props.get("h_donors", 0)
        )
    },
    {
        "name": "total_h_bonding_capacity",
        "description": "Total hydrogen bonding capacity (donors + acceptors)",
        "calculation": lambda props: props.get("h_donors", 0) + props.get("h_acceptors", 0)
    },
    {
        "name": "polarity_index",
        "description": "Ratio of TPSA to molecular surface area",
        "calculation": lambda props: props.get("tpsa", 0) / (props.get("molecular_weight", 100) ** (2/3) * 10)
    },
    {
        "name": "membrane_interaction_score",
        "description": "Estimated ability to interact with cell membranes",
        "calculation": lambda props: (
            (5.0 - abs(props.get("logp", 0) - 1.5)) / 5.0 * 
            (1.0 if props.get("molecular_weight", 0) < 400 else 0.7) *
            min(1.0, props.get("h_donors", 0) / 3.0 + props.get("h_acceptors", 0) / 5.0)
        )
    },
    {
        "name": "ice_interaction_potential",
        "description": "Estimated ability to disrupt ice crystal formation",
        "calculation": lambda props: (
            min(1.5, props.get("h_donors", 0) / 2.0) * 
            (0.5 + min(0.5, props.get("tpsa", 0) / 150.0)) *
            (1.0 if props.get("molecular_weight", 0) < 500 else 
             0.7 if props.get("molecular_weight", 0) < 1000 else 0.3)
        )
    },
    {
        "name": "vitrification_potential",
        "description": "Estimated ability to promote vitrification",
        "calculation": lambda props: (
            min(1.0, props.get("h_donors", 0) / 3.0) *
            min(1.0, props.get("h_acceptors", 0) / 6.0) *
            (1.0 if -2.0 <= props.get("logp", 0) <= 2.0 else 0.5) *
            (1.2 if props.get("ring_count", 0) > 0 else 0.8)
        )
    },
    {
        "name": "estimated_toxicity",
        "description": "Simplified toxicity estimation",
        "calculation": lambda props: (
            0.5 + 
            (0.25 if props.get("logp", 0) > 4.0 else 0.0) +
            (0.25 if props.get("molecular_weight", 0) > 500 else 0.0) -
            min(0.25, (props.get("h_donors", 0) + props.get("h_acceptors", 0)) / 20.0)
        )
    },
]

def connect_to_supabase() -> Client:
    """Connect to the Supabase database."""
    # Get URL and key from environment variables
    supabase_url = os.environ.get("SUPABASE_URL")
    service_role_key = os.environ.get("SUPABASE_SERVICE_ROLE_KEY")
    
    # If not in environment, try to load from .env file
    if not supabase_url or not service_role_key:
        try:
            load_dotenv()
            supabase_url = os.environ.get("SUPABASE_URL")
            service_role_key = os.environ.get("SUPABASE_SERVICE_ROLE_KEY")
        except Exception:
            logger.warning("Could not load .env file")
    
    # If still not available, prompt user
    if not supabase_url:
        supabase_url = input("Enter Supabase URL (https://<project>.supabase.co): ")
    
    if not service_role_key:
        import getpass
        service_role_key = getpass.getpass("Enter Supabase service role key: ")
    
    if not supabase_url or not service_role_key:
        logger.error("Missing Supabase credentials.")
        sys.exit(1)
    
    # Save to environment for future use
    os.environ["SUPABASE_URL"] = supabase_url
    os.environ["SUPABASE_SERVICE_ROLE_KEY"] = service_role_key
    
    try:
        # Use service role key for admin privileges to bypass RLS
        client = create_client(supabase_url, service_role_key)
        logger.info(f"Successfully connected to Supabase at {supabase_url}")
        return client
    except Exception as e:
        logger.error(f"Failed to connect to Supabase: {e}")
        sys.exit(1)

def get_calculation_method_id(client: Client) -> str:
    """Get or create the RDKit calculation method ID."""
    # Check if the RDKit method already exists
    rdkit_status = get_rdkit_status()
    method_name = f"RDKit-Wrapper_{rdkit_status['rdkit_version'] or 'Mock'}"
    
    response = client.table("calculation_methods").select("id").eq("name", method_name).execute()
    
    if response.data and len(response.data) > 0:
        return response.data[0]["id"]
    
    # Create a new method if it doesn't exist
    method_id = str(uuid.uuid4())
    response = client.table("calculation_methods").insert({
        "id": method_id,
        "name": method_name,
        "description": f"Molecular properties calculated using the RDKit wrapper {'with authentic RDKit' if RDKIT_AVAILABLE else 'with mock implementation'}",
        "version": rdkit_status['rdkit_version'] or "Mock",
        "method_type": "computational",
        "reference": "https://www.rdkit.org/",
        "parameters": {
            "version": rdkit_status['rdkit_version'] or "Mock",
            "calculation_type": "predictive"
        }
    }).execute()
    
    if response.data and len(response.data) > 0:
        return response.data[0]["id"]
    else:
        logger.error(f"Failed to create calculation method: {response.error}")
        return None

def get_property_types(client: Client) -> dict:
    """Get all property types needed for RDKit calculations."""
    # Get existing property types
    response = client.table("property_types").select("id, name").execute()
    
    if response.error:
        logger.error(f"Failed to get property types: {response.error}")
        return {}
    
    # Map property names to their IDs
    property_map = {item["name"]: item["id"] for item in response.data}
    
    # Standard RDKit properties
    standard_properties = [
        "molecular_weight", "logp", "tpsa", "h_donors", "h_acceptors",
        "rotatable_bonds", "ring_count", "aromatic_ring_count", "heavy_atom_count",
        "fraction_csp3", "num_stereocenters"
    ]
    
    # Cryoprotectant-specific properties
    cryo_properties = [prop["name"] for prop in CRYOPROTECTANT_PROPERTIES]
    cryo_properties.append("cryoprotectant_score")  # Add the final score
    
    # Ensure all property types exist
    for prop_name in standard_properties + cryo_properties:
        if prop_name not in property_map:
            logger.warning(f"Property type '{prop_name}' not found in database.")
    
    return property_map

def get_molecules(client: Client, args: argparse.Namespace) -> list:
    """Get molecules from the database based on command-line args."""
    query = client.table("molecules").select("id, name, smiles")
    
    # Apply filters based on args
    if args.known:
        # Known cryoprotectants search pattern
        known_terms = [
            'glycerol', 'dmso', 'dimethyl sulfoxide', 'ethylene glycol', 
            'propylene glycol', 'trehalose', 'sucrose', 'mannitol', 
            'dextran', 'ficoll', 'polyvinylpyrrolidone', 'hydroxyethyl starch',
            'formamide', 'methanol', 'glucose', 'sorbitol'
        ]
        
        # Start with an empty filter
        filter_query = client.table("molecules").select("id, name, smiles")
        
        # Add each search term with OR logic
        for i, term in enumerate(known_terms):
            if i == 0:
                filter_query = filter_query.ilike("name", f"%{term}%")
            else:
                filter_query = filter_query.or_(f"name.ilike.%{term}%")
        
        # Execute the query
        response = filter_query.execute()
    
    elif args.sample:
        # Random sample of molecules
        query = query.limit(int(args.sample)).order("random()")
        response = query.execute()
    
    elif args.limit:
        # Limit to a specific number of molecules
        query = query.limit(int(args.limit))
        response = query.execute()
    
    else:
        # No filters, get all molecules
        response = query.execute()
    
    if response.error:
        logger.error(f"Failed to get molecules: {response.error}")
        return []
    
    return response.data

def calculate_cryoprotectant_score(all_props: dict) -> float:
    """Calculate the overall cryoprotectant score based on all properties."""
    return (
        # Weighted combination of key factors
        all_props.get("membrane_interaction_score", 0) * 0.25 +
        all_props.get("ice_interaction_potential", 0) * 0.25 +
        all_props.get("vitrification_potential", 0) * 0.25 + 
        # Inverse of toxicity (1 - toxicity) with lower weight
        (1.0 - all_props.get("estimated_toxicity", 0.5)) * 0.15 +
        # Bonus for balanced hydrogen bonding
        (0.1 if 0.5 <= all_props.get("h_bond_donor_acceptor_ratio", 0) <= 2.0 else 0.0)
    ) * 10  # Scale to 0-10 range

def store_property(client: Client, molecule_id: str, property_type_id: str, 
                  property_name: str, numeric_value: float, 
                  calculation_method_id: str) -> bool:
    """Store a property value in the database."""
    if numeric_value is None:
        return False
    
    try:
        # Check if property already exists
        response = client.table("molecular_properties")\
            .select("id")\
            .eq("molecule_id", molecule_id)\
            .eq("property_type_id", property_type_id)\
            .execute()
        
        if response.data and len(response.data) > 0:
            # Update existing property
            existing_id = response.data[0]["id"]
            response = client.table("molecular_properties")\
                .update({"numeric_value": numeric_value, "updated_at": datetime.now().isoformat()})\
                .eq("id", existing_id)\
                .execute()
            return True
        else:
            # Insert new property
            property_id = str(uuid.uuid4())
            response = client.table("molecular_properties").insert({
                "id": property_id,
                "molecule_id": molecule_id,
                "property_type_id": property_type_id,
                "numeric_value": numeric_value,
                "property_name": property_name,
                "property_type": "numeric",
                "source": "rdkit_wrapper_calculation"
            }).execute()
            
            return True
            
    except Exception as e:
        logger.error(f"Error storing property {property_name} for molecule {molecule_id}: {e}")
        return False

def calculate_and_store_properties(client: Client, molecule: dict, 
                                  property_types: dict, 
                                  calculation_method_id: str) -> dict:
    """Calculate and store all RDKit properties for a molecule."""
    molecule_id = molecule["id"]
    smiles = molecule["smiles"]
    
    if not smiles:
        logger.warning(f"Skipping molecule {molecule['name']} (ID: {molecule_id}) - No SMILES available")
        return {"status": "skipped", "reason": "no_smiles"}
    
    # Create molecule from SMILES
    mol = create_molecule_from_smiles(smiles)
    if not mol:
        logger.warning(f"Failed to create molecule from SMILES for {molecule['name']} (ID: {molecule_id})")
        return {"status": "skipped", "reason": "invalid_smiles"}
    
    # Calculate standard properties
    try:
        # Get RDKit properties
        props = calculate_properties(mol)
        
        # Store results
        results = {
            "stored_properties": 0,
            "calculated_properties": {}
        }
        
        # Store standard RDKit properties
        standard_properties = [
            "molecular_weight", "logp", "tpsa", "h_donors", "h_acceptors",
            "rotatable_bonds", "ring_count", "aromatic_ring_count", "heavy_atom_count"
        ]
        
        for prop_name in standard_properties:
            if prop_name in props and props[prop_name] is not None:
                prop_type_id = property_types.get(prop_name)
                if prop_type_id:
                    if store_property(client, molecule_id, prop_type_id, 
                                      prop_name, props[prop_name], 
                                      calculation_method_id):
                        results["stored_properties"] += 1
                    
                    # Save for derived calculations
                    results["calculated_properties"][prop_name] = props[prop_name]
        
        # Calculate and store cryoprotectant-specific properties
        all_props = results["calculated_properties"].copy()
        for prop in CRYOPROTECTANT_PROPERTIES:
            try:
                prop_name = prop["name"]
                prop_value = prop["calculation"](all_props)
                
                prop_type_id = property_types.get(prop_name)
                if prop_type_id:
                    if store_property(client, molecule_id, prop_type_id, 
                                     prop_name, prop_value, 
                                     calculation_method_id):
                        results["stored_properties"] += 1
                    
                    # Add to properties for later calculations
                    all_props[prop_name] = prop_value
            except Exception as e:
                logger.error(f"Error calculating {prop_name} for molecule {molecule['name']}: {e}")
        
        # Calculate and store the overall cryoprotectant score
        cryo_score = calculate_cryoprotectant_score(all_props)
        cryo_score_type_id = property_types.get("cryoprotectant_score")
        
        if cryo_score_type_id:
            if store_property(client, molecule_id, cryo_score_type_id, 
                             "cryoprotectant_score", cryo_score, 
                             calculation_method_id):
                results["stored_properties"] += 1
                all_props["cryoprotectant_score"] = cryo_score
        
        return {
            "status": "success", 
            "results": results,
            "molecule_name": molecule['name'],
            "cryo_score": cryo_score
        }
        
    except Exception as e:
        logger.error(f"Error calculating properties for molecule {molecule['name']} (ID: {molecule_id}): {e}")
        return {"status": "error", "error": str(e)}

def log_calculation_results(client: Client, results: dict, calculation_method_id: str):
    """Log the calculation results to the database."""
    try:
        log_id = str(uuid.uuid4())
        response = client.table("scientific_data_audit").insert({
            "id": log_id,
            "operation": "RDKIT_PROPERTY_CALCULATION",
            "timestamp": datetime.now().isoformat(),
            "table_name": "molecular_properties",
            "application_context": "rdkit_property_calculator.py",
            "new_data": results
        }).execute()
        
        logger.info(f"Logged calculation results to audit table with ID: {log_id}")
    except Exception as e:
        logger.warning(f"Failed to log results to audit table: {e}")

def main():
    """Main function for RDKit property calculation."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Calculate and store RDKit molecular properties")
    parser.add_argument("--known", action="store_true", help="Process only known cryoprotectants")
    parser.add_argument("--sample", type=int, help="Process a random sample of N molecules")
    parser.add_argument("--limit", type=int, help="Limit processing to N molecules")
    parser.add_argument("--url", type=str, help="Supabase URL (overrides environment variable)")
    parser.add_argument("--key", type=str, help="Supabase service role key (overrides environment variable)")
    args = parser.parse_args()
    
    # Set credentials from command line if provided
    if args.url:
        os.environ["SUPABASE_URL"] = args.url
    if args.key:
        os.environ["SUPABASE_SERVICE_ROLE_KEY"] = args.key
    
    # Print RDKit status
    rdkit_status = get_rdkit_status()
    logger.info(f"RDKit Status: {'Available' if rdkit_status['rdkit_available'] else 'Not Available'}")
    if rdkit_status['rdkit_available']:
        logger.info(f"RDKit Version: {rdkit_status['rdkit_version']}")
    else:
        logger.warning("Using mock implementation for calculations")
    
    try:
        # Connect to Supabase
        logger.info("Connecting to Supabase database...")
        client = connect_to_supabase()
        
        # Get or create calculation method
        logger.info("Setting up calculation method...")
        calculation_method_id = get_calculation_method_id(client)
        logger.info(f"Using calculation method ID: {calculation_method_id}")
        
        # Get all property types
        logger.info("Loading property types...")
        property_types = get_property_types(client)
        logger.info(f"Found {len(property_types)} property types")
        
        # Get molecules based on arguments
        logger.info("Fetching molecules...")
        molecules = get_molecules(client, args)
        logger.info(f"Processing {len(molecules)} molecules")
        
        if not molecules:
            logger.error("No molecules found to process. Check your filters or database connection.")
            sys.exit(1)
        
        # Initialize counters
        successful = 0
        skipped = 0
        errors = 0
        high_scoring = []
        
        # Process all molecules with a progress bar
        for molecule in tqdm(molecules, desc="Processing molecules"):
            try:
                result = calculate_and_store_properties(
                    client, molecule, property_types, calculation_method_id
                )
                
                if result["status"] == "success":
                    successful += 1
                    
                    # Track high-scoring molecules
                    if result.get("cryo_score", 0) > 7.0:
                        high_scoring.append({
                            "name": molecule["name"],
                            "score": result["cryo_score"]
                        })
                        
                elif result["status"] == "skipped":
                    skipped += 1
                else:
                    errors += 1
            except Exception as e:
                logger.error(f"Error processing molecule {molecule.get('name', 'unknown')}: {e}")
                errors += 1
        
        # Log summary
        summary = {
            "total_molecules": len(molecules),
            "successful": successful,
            "skipped": skipped,
            "errors": errors,
            "high_scoring_count": len(high_scoring),
            "rdkit_available": rdkit_status['rdkit_available'],
            "rdkit_version": rdkit_status['rdkit_version'],
            "calculation_method_id": calculation_method_id,
            "high_scoring_molecules": high_scoring[:10]  # Show top 10
        }
        
        # Log results to database
        try:
            log_calculation_results(client, summary, calculation_method_id)
        except Exception as e:
            logger.error(f"Failed to log calculation results: {e}")
        
        logger.info(f"\n{'-'*40}")
        logger.info("Calculation Summary:")
        logger.info(f"Total molecules processed: {len(molecules)}")
        logger.info(f"Successful calculations: {successful}")
        logger.info(f"Skipped molecules: {skipped}")
        logger.info(f"Failed calculations: {errors}")
        
        if high_scoring:
            logger.info(f"\nTop Cryoprotectant Candidates:")
            for i, mol in enumerate(sorted(high_scoring, key=lambda x: x["score"], reverse=True)[:5]):
                logger.info(f"{i+1}. {mol['name']} (Score: {mol['score']:.2f})")
        
        logger.info(f"{'-'*40}")
        
    except KeyboardInterrupt:
        logger.warning("Process interrupted by user.")
        sys.exit(1)
    except supabase.Client.AuthError as e:
        logger.error(f"Authentication error: {e}")
        logger.error("This is likely due to an invalid service role key or URL.")
        logger.error("Please check your credentials in the .env file or use --url and --key parameters.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error processing molecules: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()