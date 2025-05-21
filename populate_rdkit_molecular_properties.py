#!/usr/bin/env python3
"""
Populate molecular properties using RDKit calculations

This script calculates and populates molecular properties for all molecules
in the database using the new rdkit_wrapper. The wrapper provides a consistent
interface to RDKit with fallback to simpler calculations when RDKit isn't available.

The script focuses on properties relevant to cryoprotectant viability prediction:
- Basic molecular properties (weight, LogP, etc.)
- Hydrogen bonding capabilities
- Molecular polarity metrics
- Structural features that influence cryoprotection

The calculated properties are stored in the molecular_properties table.
"""

import os
import sys
import uuid
import logging
import psycopg2
import argparse
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
from tqdm import tqdm
import json

# Import our RDKit wrapper
from rdkit_wrapper import (
    create_molecule_from_smiles,
    calculate_properties,
    get_rdkit_status,
    generate_fingerprint,
    RDKIT_AVAILABLE
)

# Configure logging
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

# Cryoprotectant-specific properties to calculate beyond standard RDKit properties
CRYOPROTECTANT_PROPERTIES = [
    # Hydrogen bonding capacity is critical for cryoprotection
    {
        "name": "h_bond_donor_acceptor_ratio",
        "description": "Ratio of H-bond donors to acceptors, indicator of hydrogen bonding symmetry",
        "data_type": "numeric",
        "units": None,
        "calculation": lambda props: (
            props.get("h_donors", 0) / props.get("h_acceptors", 1)
            if props.get("h_acceptors", 0) > 0 else props.get("h_donors", 0)
        )
    },
    {
        "name": "total_h_bonding_capacity",
        "description": "Total hydrogen bonding capacity (donors + acceptors)",
        "data_type": "numeric",
        "units": None,
        "calculation": lambda props: props.get("h_donors", 0) + props.get("h_acceptors", 0)
    },
    # Polarity and membrane interaction metrics
    {
        "name": "polarity_index",
        "description": "Ratio of TPSA to molecular surface area, indicator of overall polarity",
        "data_type": "numeric",
        "units": None,
        "calculation": lambda props: props.get("tpsa", 0) / (props.get("molecular_weight", 100) ** (2/3) * 10)
        # Using molecular weight^(2/3) as approximate surface area scaling
    },
    {
        "name": "membrane_interaction_score",
        "description": "Estimated ability to interact with cell membranes based on LogP and structure",
        "data_type": "numeric", 
        "units": None,
        "calculation": lambda props: (
            # Balanced LogP (not too hydrophobic or hydrophilic) is ideal for membrane interaction
            (5.0 - abs(props.get("logp", 0) - 1.5)) / 5.0 * 
            # Small molecules with balanced polarity work better
            (1.0 if props.get("molecular_weight", 0) < 400 else 0.7) *
            # Hydrogen bonding capacity helps stabilize membranes
            min(1.0, props.get("h_donors", 0) / 3.0 + props.get("h_acceptors", 0) / 5.0)
        )
    },
    # Ice formation inhibition capacity
    {
        "name": "ice_interaction_potential",
        "description": "Estimated ability to interact with and disrupt ice crystal formation",
        "data_type": "numeric",
        "units": None, 
        "calculation": lambda props: (
            # Hydroxyl groups are excellent ice disruptors
            min(1.5, props.get("h_donors", 0) / 2.0) * 
            # TPSA represents polar surface that can interact with water
            (0.5 + min(0.5, props.get("tpsa", 0) / 150.0)) *
            # Reasonable molecular size needed (not too large)
            (1.0 if props.get("molecular_weight", 0) < 500 else 
             0.7 if props.get("molecular_weight", 0) < 1000 else 0.3)
        )
    },
    # Vitrification potential
    {
        "name": "vitrification_potential", 
        "description": "Estimated ability to promote vitrification (glass formation) during freezing",
        "data_type": "numeric",
        "units": None,
        "calculation": lambda props: (
            # Complex calculation based on several factors that contribute to vitrification
            # 1. Hydroxyl groups help with water interaction
            min(1.0, props.get("h_donors", 0) / 3.0) *
            # 2. Larger molecules with many H-bond acceptors promote networking
            min(1.0, props.get("h_acceptors", 0) / 6.0) *
            # 3. Moderate LogP (not too hydrophobic)
            (1.0 if -2.0 <= props.get("logp", 0) <= 2.0 else 0.5) *
            # 4. Ring structures can inhibit crystallization
            (1.2 if props.get("ring_count", 0) > 0 else 0.8)
        )
    },
    # Toxicity estimation - simplified approach
    {
        "name": "estimated_toxicity",
        "description": "Simplified toxicity estimation based on molecular properties",
        "data_type": "numeric",
        "units": None,
        "calculation": lambda props: (
            # Start with moderate score
            0.5 + 
            # Very high LogP indicates higher toxicity
            (0.25 if props.get("logp", 0) > 4.0 else 0.0) +
            # Very high molecular weight can indicate reduced clearance
            (0.25 if props.get("molecular_weight", 0) > 500 else 0.0) -
            # More H-bonding generally means better water solubility
            min(0.25, (props.get("h_donors", 0) + props.get("h_acceptors", 0)) / 20.0)
        )
    },
    # Overall cryoprotectant score
    {
        "name": "cryoprotectant_score",
        "description": "Overall estimated effectiveness as a cryoprotectant",
        "data_type": "numeric",
        "units": None,
        "calculated_at_end": True,
        "calculation": lambda all_props: (
            # Weighted combination of key factors
            all_props.get("membrane_interaction_score", 0) * 0.25 +
            all_props.get("ice_interaction_potential", 0) * 0.25 +
            all_props.get("vitrification_potential", 0) * 0.25 + 
            # Inverse of toxicity (1 - toxicity) with lower weight
            (1.0 - all_props.get("estimated_toxicity", 0.5)) * 0.15 +
            # Bonus for balanced hydrogen bonding
            (0.1 if 0.5 <= all_props.get("h_bond_donor_acceptor_ratio", 0) <= 2.0 else 0.0)
        ) * 10  # Scale to 0-10 range
    }
]

# Standard RDKit properties to store in the database
STANDARD_RDKIT_PROPERTIES = [
    {"name": "molecular_weight", "data_type": "numeric", "units": "g/mol"},
    {"name": "exact_mass", "data_type": "numeric", "units": "g/mol"},
    {"name": "logp", "data_type": "numeric", "units": None},
    {"name": "tpsa", "data_type": "numeric", "units": "Å²"},
    {"name": "h_donors", "data_type": "numeric", "units": None},
    {"name": "h_acceptors", "data_type": "numeric", "units": None},
    {"name": "rotatable_bonds", "data_type": "numeric", "units": None},
    {"name": "ring_count", "data_type": "numeric", "units": None},
    {"name": "aromatic_ring_count", "data_type": "numeric", "units": None},
    {"name": "heavy_atom_count", "data_type": "numeric", "units": None},
    {"name": "fraction_csp3", "data_type": "numeric", "units": None},
    {"name": "num_stereocenters", "data_type": "numeric", "units": None}
]

def connect_to_db():
    """Connect to the database using environment variables."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    
    # Ensure required parameters are present
    if not all([db_params['host'], db_params['user'], db_params['password']]):
        logger.error("Missing required database connection parameters in environment variables.")
        logger.error("Make sure you have a valid .env file with SUPABASE_DB_* variables.")
        sys.exit(1)
    
    try:
        conn = psycopg2.connect(**db_params)
        return conn
    except psycopg2.Error as e:
        logger.error(f"Database connection error: {e}")
        sys.exit(1)

def get_calculation_method_id(conn):
    """Get or create the calculation method for RDKit-based calculations."""
    method_name = f"RDKit-Wrapper_{get_rdkit_status()['rdkit_version'] or 'Mock'}"
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Check if method exists
        cursor.execute("""
            SELECT id FROM calculation_methods WHERE name = %s
        """, (method_name,))
        
        result = cursor.fetchone()
        if result:
            return result['id']
        
        # Create new method
        method_id = str(uuid.uuid4())
        cursor.execute("""
            INSERT INTO calculation_methods (
                id, name, description, version, created_at, updated_at
            ) VALUES (
                %s, %s, %s, %s, NOW(), NOW()
            ) RETURNING id
        """, (
            method_id,
            method_name,
            f"Molecular properties calculated using the RDKit wrapper {'with authentic RDKit' if RDKIT_AVAILABLE else 'with mock implementation'}",
            get_rdkit_status()['rdkit_version'] or 'Mock'
        ))
        
        return cursor.fetchone()['id']

def get_or_create_property_type(conn, prop_name, prop_type, description, units):
    """Get or create a property type in the database."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Check if property type exists
        cursor.execute("""
            SELECT id FROM property_types WHERE name = %s
        """, (prop_name,))
        
        result = cursor.fetchone()
        if result:
            return result['id']
        
        # Create new property type
        prop_id = str(uuid.uuid4())
        cursor.execute("""
            INSERT INTO property_types (
                id, name, description, units, data_type, created_at, updated_at
            ) VALUES (
                %s, %s, %s, %s, %s, NOW(), NOW()
            ) RETURNING id
        """, (
            prop_id,
            prop_name,
            description,
            units,
            prop_type
        ))
        
        return cursor.fetchone()['id']

def get_property_types(conn, property_list):
    """Get or create all property types needed for the calculations."""
    property_type_ids = {}
    
    for prop in property_list:
        prop_id = get_or_create_property_type(
            conn, 
            prop['name'], 
            prop['data_type'], 
            prop.get('description', prop['name']), 
            prop.get('units', None)
        )
        property_type_ids[prop['name']] = prop_id
    
    return property_type_ids

def get_molecules(conn, batch_size=100, where_clause=""):
    """Get all molecules from the database in batches."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Count total molecules for progress tracking
        cursor.execute(f"SELECT COUNT(*) FROM molecules {where_clause}")
        total = cursor.fetchone()['count']
        
        # Initialize offset
        offset = 0
        
        # Create progress bar
        with tqdm(total=total, desc="Processing molecules") as pbar:
            while True:
                cursor.execute(f"""
                    SELECT id, name, smiles, pubchem_cid, molecular_formula 
                    FROM molecules 
                    {where_clause}
                    ORDER BY id
                    LIMIT %s OFFSET %s
                """, (batch_size, offset))
                
                batch = cursor.fetchall()
                if not batch:
                    break
                
                yield batch
                
                offset += batch_size
                pbar.update(len(batch))

def store_property(conn, molecule_id, property_type_id, value, property_name, calculation_method_id):
    """Store a calculated property in the database."""
    if value is None:
        return False
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        try:
            # Check if property already exists
            cursor.execute("""
                SELECT id FROM molecular_properties 
                WHERE molecule_id = %s AND property_type_id = %s
            """, (molecule_id, property_type_id))
            
            existing = cursor.fetchone()
            
            if existing:
                # Update existing property
                cursor.execute("""
                    UPDATE molecular_properties SET
                        numeric_value = %s,
                        updated_at = NOW()
                    WHERE id = %s
                """, (value, existing['id']))
                return "updated"
            else:
                # Insert new property
                prop_id = str(uuid.uuid4())
                cursor.execute("""
                    INSERT INTO molecular_properties (
                        id, molecule_id, property_type_id, numeric_value, 
                        created_at, updated_at
                    ) VALUES (
                        %s, %s, %s, %s, NOW(), NOW()
                    )
                """, (
                    prop_id,
                    molecule_id,
                    property_type_id,
                    value
                ))
                
                # Also store in predictions table if method_id is provided
                if calculation_method_id:
                    pred_id = str(uuid.uuid4())
                    cursor.execute("""
                        INSERT INTO predictions (
                            id, molecule_id, property_type_id, 
                            calculation_method_id, numeric_value,
                            created_at, updated_at
                        ) VALUES (
                            %s, %s, %s, %s, %s, NOW(), NOW()
                        ) ON CONFLICT (molecule_id, property_type_id, calculation_method_id) 
                        DO UPDATE SET 
                            numeric_value = EXCLUDED.numeric_value,
                            updated_at = NOW()
                    """, (
                        pred_id,
                        molecule_id,
                        property_type_id,
                        calculation_method_id,
                        value
                    ))
                
                return "inserted"
            
        except Exception as e:
            logger.error(f"Error storing property {property_name} for molecule {molecule_id}: {e}")
            return False

def calculate_and_store_properties(conn, molecule, property_type_ids, calculation_method_id):
    """Calculate and store all properties for a molecule."""
    molecule_id = molecule['id']
    smiles = molecule['smiles']
    
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
        props = calculate_properties(mol)
        
        # Store results
        results = {
            "inserted": 0,
            "updated": 0,
            "failed": 0,
            "calculated_properties": {}
        }
        
        # Process standard RDKit properties
        for prop in STANDARD_RDKIT_PROPERTIES:
            prop_name = prop['name']
            if prop_name in props and props[prop_name] is not None:
                prop_type_id = property_type_ids.get(prop_name)
                if prop_type_id:
                    result = store_property(
                        conn, 
                        molecule_id, 
                        prop_type_id, 
                        props[prop_name], 
                        prop_name,
                        calculation_method_id
                    )
                    
                    if result == "inserted":
                        results["inserted"] += 1
                    elif result == "updated":
                        results["updated"] += 1
                    else:
                        results["failed"] += 1
                    
                    # Save for use in derived calculations
                    results["calculated_properties"][prop_name] = props[prop_name]
        
        # Calculate and store cryoprotectant-specific properties
        all_props = results["calculated_properties"].copy()
        for prop in CRYOPROTECTANT_PROPERTIES:
            # Skip properties that need to be calculated at the end
            if prop.get('calculated_at_end', False):
                continue
                
            try:
                prop_name = prop['name']
                prop_value = prop['calculation'](all_props)
                
                prop_type_id = property_type_ids.get(prop_name)
                if prop_type_id:
                    result = store_property(
                        conn, 
                        molecule_id, 
                        prop_type_id, 
                        prop_value, 
                        prop_name,
                        calculation_method_id
                    )
                    
                    if result == "inserted":
                        results["inserted"] += 1
                    elif result == "updated":
                        results["updated"] += 1
                    else:
                        results["failed"] += 1
                    
                    # Add to properties for later calculations
                    all_props[prop_name] = prop_value
            except Exception as e:
                logger.error(f"Error calculating {prop_name} for molecule {molecule['name']}: {e}")
                results["failed"] += 1
        
        # Calculate properties that depend on other calculated properties
        for prop in CRYOPROTECTANT_PROPERTIES:
            if prop.get('calculated_at_end', False):
                try:
                    prop_name = prop['name']
                    prop_value = prop['calculation'](all_props)
                    
                    prop_type_id = property_type_ids.get(prop_name)
                    if prop_type_id:
                        result = store_property(
                            conn, 
                            molecule_id, 
                            prop_type_id, 
                            prop_value, 
                            prop_name,
                            calculation_method_id
                        )
                        
                        if result == "inserted":
                            results["inserted"] += 1
                        elif result == "updated":
                            results["updated"] += 1
                        else:
                            results["failed"] += 1
                except Exception as e:
                    logger.error(f"Error calculating final {prop_name} for molecule {molecule['name']}: {e}")
                    results["failed"] += 1
        
        return {
            "status": "success", 
            "results": results,
            "molecule_name": molecule['name']
        }
        
    except Exception as e:
        logger.error(f"Error calculating properties for molecule {molecule['name']} (ID: {molecule_id}): {e}")
        return {"status": "error", "error": str(e)}

def log_results_to_database(conn, summary):
    """Log calculation results to an audit table for tracking."""
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
                INSERT INTO scientific_data_audit (
                    operation,
                    timestamp,
                    table_name,
                    user_id,
                    old_data,
                    new_data,
                    application_context
                ) VALUES (
                    'RDKIT_PROPERTY_CALCULATION',
                    NOW(),
                    'molecular_properties',
                    NULL,
                    NULL,
                    %s,
                    'populate_rdkit_molecular_properties.py'
                )
            """, (json.dumps(summary),))
    except Exception as e:
        logger.warning(f"Failed to log results to audit table: {e}")

def main():
    """Main function for RDKit property calculation."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Calculate and store RDKit molecular properties")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for processing molecules")
    parser.add_argument("--where", type=str, default="", help="Optional WHERE clause for molecule selection")
    parser.add_argument("--commit-interval", type=int, default=10, help="Number of batches to process before committing")
    args = parser.parse_args()
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False  # Start transaction mode
    
    # Print RDKit status
    rdkit_status = get_rdkit_status()
    logger.info(f"RDKit Status: {'Available' if rdkit_status['rdkit_available'] else 'Not Available'}")
    if rdkit_status['rdkit_available']:
        logger.info(f"RDKit Version: {rdkit_status['rdkit_version']}")
    else:
        logger.warning("Using mock implementation for calculations")
    
    try:
        # Get or create calculation method
        calculation_method_id = get_calculation_method_id(conn)
        logger.info(f"Using calculation method ID: {calculation_method_id}")
        
        # Get all property types
        all_properties = STANDARD_RDKIT_PROPERTIES + CRYOPROTECTANT_PROPERTIES
        property_type_ids = get_property_types(conn, all_properties)
        logger.info(f"Using {len(property_type_ids)} property types")
        
        # Process molecules in batches
        where_clause = f"WHERE {args.where}" if args.where else ""
        
        # Initialize counters
        total_molecules = 0
        successful = 0
        skipped = 0
        errors = 0
        properties_inserted = 0
        properties_updated = 0
        batch_count = 0
        
        # Process all molecules
        for batch in get_molecules(conn, args.batch_size, where_clause):
            batch_results = []
            
            for molecule in batch:
                result = calculate_and_store_properties(
                    conn, molecule, property_type_ids, calculation_method_id
                )
                
                batch_results.append(result)
                total_molecules += 1
                
                if result["status"] == "success":
                    successful += 1
                    properties_inserted += result["results"]["inserted"]
                    properties_updated += result["results"]["updated"]
                elif result["status"] == "skipped":
                    skipped += 1
                else:
                    errors += 1
            
            # Commit periodically to avoid long transactions
            batch_count += 1
            if batch_count % args.commit_interval == 0:
                conn.commit()
                logger.info(f"Committed after {batch_count} batches. Processed {total_molecules} molecules so far.")
        
        # Commit final changes
        conn.commit()
        
        # Log summary
        summary = {
            "total_molecules": total_molecules,
            "successful": successful,
            "skipped": skipped,
            "errors": errors,
            "properties_inserted": properties_inserted,
            "properties_updated": properties_updated,
            "rdkit_available": rdkit_status['rdkit_available'],
            "rdkit_version": rdkit_status['rdkit_version'],
            "calculation_method_id": calculation_method_id
        }
        
        # Log results to database
        log_results_to_database(conn, summary)
        
        logger.info(f"\n{'-'*40}")
        logger.info("Calculation Summary:")
        logger.info(f"Total molecules processed: {total_molecules}")
        logger.info(f"Successful calculations: {successful}")
        logger.info(f"Skipped molecules: {skipped}")
        logger.info(f"Failed calculations: {errors}")
        logger.info(f"Properties inserted: {properties_inserted}")
        logger.info(f"Properties updated: {properties_updated}")
        logger.info(f"{'-'*40}")
        
    except KeyboardInterrupt:
        logger.warning("Process interrupted by user. Rolling back changes.")
        conn.rollback()
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error processing molecules: {e}")
        conn.rollback()
        sys.exit(1)
    finally:
        conn.close()

if __name__ == "__main__":
    main()