#!/usr/bin/env python3
"""
Database Update Scheduled Job

This script is designed to run on a schedule to update database data,
including molecular properties, statistics, and cryoprotection scores.
"""

import os
import sys
import json
import logging
from datetime import datetime
import traceback

# Add project directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'logs',
            f'db_update_{datetime.now().strftime("%Y%m%d")}.log'
        )),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('database_update')

# Import database utilities
try:
    from database.adapter import get_database_connection
    from api.property_explorer_resources import get_property_statistics
    from api.utils import get_supabase_client
except ImportError as e:
    logger.error(f"Error importing required modules: {e}")
    sys.exit(1)

def ensure_log_directory():
    """Ensure log directory exists"""
    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

def update_property_statistics():
    """Update property statistics"""
    logger.info("Updating property statistics")
    try:
        supabase = get_supabase_client()
        
        # Get all property types
        property_types_response = supabase.from_("property_types").select("id, name").execute()
        if not hasattr(property_types_response, 'data'):
            logger.error("Failed to retrieve property types")
            return False
        
        property_types = property_types_response.data
        
        # Calculate and update statistics for each property type
        updated_count = 0
        for prop_type in property_types:
            try:
                # Calculate statistics
                stats = get_property_statistics(supabase, prop_type['id'])
                
                # Update or create statistics record
                stats_response = supabase.from_("property_statistics").upsert({
                    "property_type_id": prop_type['id'],
                    "min_value": stats.get('min'),
                    "max_value": stats.get('max'),
                    "mean_value": stats.get('mean'),
                    "median_value": stats.get('median'),
                    "std_dev": stats.get('std_dev'),
                    "count": stats.get('count'),
                    "distribution": json.dumps(stats.get('distribution', [])),
                    "last_updated": datetime.now().isoformat()
                }).execute()
                
                updated_count += 1
                logger.info(f"Updated statistics for property {prop_type['name']}")
            except Exception as e:
                logger.error(f"Error updating statistics for property {prop_type['name']}: {e}")
                continue
        
        logger.info(f"Updated statistics for {updated_count} properties")
        return True
    except Exception as e:
        logger.error(f"Error updating property statistics: {e}")
        logger.error(traceback.format_exc())
        return False

def update_cryoprotection_scores():
    """Update cryoprotection scores for molecules"""
    logger.info("Updating cryoprotection scores")
    try:
        supabase = get_supabase_client()
        
        # Get molecules without cryoprotection scores
        query = """
        SELECT 
            m.id 
        FROM 
            molecules m
        LEFT JOIN 
            molecular_properties mp ON m.id = mp.molecule_id AND mp.property_type_id = (
                SELECT id FROM property_types WHERE name = 'cryoprotection_score'
            )
        WHERE 
            mp.id IS NULL
        LIMIT 100
        """
        
        response = supabase.rpc("exec_sql", {"sql_query": query}).execute()
        
        if not hasattr(response, 'data'):
            logger.error("Failed to retrieve molecules without cryoprotection scores")
            return False
        
        molecules = response.data
        
        if not molecules:
            logger.info("No molecules without cryoprotection scores found")
            return True
        
        # Calculate cryoprotection scores for each molecule
        updated_count = 0
        for molecule in molecules:
            try:
                # In a real implementation, this would use a more sophisticated
                # scoring algorithm. For this demo, we'll use a simple approach.
                
                # Get molecular properties like LogP, TPSA, etc.
                props_query = """
                SELECT 
                    pt.name, mp.value
                FROM 
                    molecular_properties mp
                JOIN 
                    property_types pt ON mp.property_type_id = pt.id
                WHERE 
                    mp.molecule_id = $1
                    AND pt.name IN ('logp', 'tpsa', 'hba', 'hbd', 'mw')
                """
                
                props_response = supabase.rpc(
                    "exec_sql", 
                    {"sql_query": props_query, "params": [molecule['id']]}
                ).execute()
                
                if not hasattr(props_response, 'data') or not props_response.data:
                    logger.warning(f"No properties found for molecule {molecule['id']}")
                    continue
                
                # Convert properties to dict
                properties = {p['name']: p['value'] for p in props_response.data}
                
                # Simple scoring algorithm (replace with actual algorithm in production)
                score = calculate_cryoprotection_score(properties)
                
                # Update or create cryoprotection score
                cryo_score_type_response = supabase.from_("property_types").select("id").eq("name", "cryoprotection_score").execute()
                
                if not hasattr(cryo_score_type_response, 'data') or not cryo_score_type_response.data:
                    logger.error("Cryoprotection score property type not found")
                    continue
                
                cryo_score_type_id = cryo_score_type_response.data[0]['id']
                
                # Insert score
                score_response = supabase.from_("molecular_properties").upsert({
                    "molecule_id": molecule['id'],
                    "property_type_id": cryo_score_type_id,
                    "value": score,
                    "source": "calculated",
                    "confidence": 0.8,
                    "last_updated": datetime.now().isoformat()
                }).execute()
                
                updated_count += 1
                logger.info(f"Updated cryoprotection score for molecule {molecule['id']}")
            except Exception as e:
                logger.error(f"Error updating cryoprotection score for molecule {molecule['id']}: {e}")
                continue
        
        logger.info(f"Updated cryoprotection scores for {updated_count} molecules")
        return True
    except Exception as e:
        logger.error(f"Error updating cryoprotection scores: {e}")
        logger.error(traceback.format_exc())
        return False

def calculate_cryoprotection_score(properties):
    """
    Calculate cryoprotection score based on molecular properties
    
    This is a simplified calculation that would be replaced with a more
    sophisticated algorithm in a real implementation.
    """
    # Default score
    score = 0.5
    
    # LogP: lower is better for cryoprotectants (more hydrophilic)
    if 'logp' in properties:
        logp = float(properties['logp'])
        # Ideal range: -2 to 1
        if -2 <= logp <= 1:
            score += 0.1
        # Penalty for very hydrophobic molecules
        elif logp > 3:
            score -= 0.1
    
    # TPSA: higher is better for cryoprotectants (more hydrogen bonding)
    if 'tpsa' in properties:
        tpsa = float(properties['tpsa'])
        # Ideal range: 40-120
        if 40 <= tpsa <= 120:
            score += 0.1
        # Bonus for very polar molecules
        elif tpsa > 120:
            score += 0.05
    
    # H-bond acceptors: more is better (up to a point)
    if 'hba' in properties:
        hba = int(properties['hba'])
        # Ideal range: 3-7
        if 3 <= hba <= 7:
            score += 0.1
        # Bonus for more acceptors
        elif hba > 7:
            score += 0.05
    
    # H-bond donors: more is better (up to a point)
    if 'hbd' in properties:
        hbd = int(properties['hbd'])
        # Ideal range: 2-5
        if 2 <= hbd <= 5:
            score += 0.1
        # Bonus for more donors
        elif hbd > 5:
            score += 0.05
    
    # Molecular weight: lower is generally better
    if 'mw' in properties:
        mw = float(properties['mw'])
        # Ideal range: 60-200
        if 60 <= mw <= 200:
            score += 0.1
        # Penalty for very large molecules
        elif mw > 400:
            score -= 0.1
    
    # Ensure score is in range 0-1
    score = max(0, min(1, score))
    
    return score

def update_molecule_relationships():
    """Update molecule relationship data"""
    logger.info("Updating molecule relationships")
    try:
        supabase = get_supabase_client()
        
        # Identify molecules that need relationship updates
        # In this demo, we'll focus on updating molecule similarities
        
        query = """
        SELECT 
            m1.id as mol1_id, 
            m2.id as mol2_id
        FROM 
            molecules m1
        CROSS JOIN 
            molecules m2
        LEFT JOIN 
            molecule_relationships mr ON 
                (mr.molecule1_id = m1.id AND mr.molecule2_id = m2.id) OR
                (mr.molecule1_id = m2.id AND mr.molecule2_id = m1.id)
        WHERE 
            m1.id < m2.id AND
            mr.id IS NULL
        LIMIT 100
        """
        
        response = supabase.rpc("exec_sql", {"sql_query": query}).execute()
        
        if not hasattr(response, 'data'):
            logger.error("Failed to retrieve molecules for relationship updates")
            return False
        
        molecule_pairs = response.data
        
        if not molecule_pairs:
            logger.info("No molecule relationships to update")
            return True
        
        # Update relationships for each molecule pair
        updated_count = 0
        for pair in molecule_pairs:
            try:
                # In a real implementation, this would calculate actual similarity
                # For this demo, we'll use a random similarity value
                import random
                similarity = random.uniform(0, 1)
                
                # Insert relationship
                rel_response = supabase.from_("molecule_relationships").insert({
                    "molecule1_id": pair['mol1_id'],
                    "molecule2_id": pair['mol2_id'],
                    "relationship_type": "similarity",
                    "relationship_value": similarity,
                    "last_updated": datetime.now().isoformat()
                }).execute()
                
                updated_count += 1
                logger.info(f"Added relationship between molecules {pair['mol1_id']} and {pair['mol2_id']}")
            except Exception as e:
                logger.error(f"Error updating relationship for molecules {pair['mol1_id']} and {pair['mol2_id']}: {e}")
                continue
        
        logger.info(f"Added {updated_count} molecule relationships")
        return True
    except Exception as e:
        logger.error(f"Error updating molecule relationships: {e}")
        logger.error(traceback.format_exc())
        return False

def perform_database_update():
    """Perform all database update tasks"""
    logger.info("Starting scheduled database update")
    
    # Run update tasks
    try:
        # Update property statistics
        stats_success = update_property_statistics()
        
        # Update cryoprotection scores
        scores_success = update_cryoprotection_scores()
        
        # Update molecule relationships
        relationships_success = update_molecule_relationships()
        
        # Log overall results
        if stats_success and scores_success and relationships_success:
            logger.info("Database update completed successfully")
            return True
        else:
            logger.warning("Database update completed with some errors")
            return False
    except Exception as e:
        logger.error(f"Error during database update: {e}")
        logger.error(traceback.format_exc())
        return False

def get_property_statistics(supabase, property_id):
    """
    Calculate statistics for a property
    
    This function is normally imported from api.property_explorer_resources,
    but we're implementing it here for this demo.
    """
    try:
        # Get property values
        response = supabase.from_("molecular_properties").select("value").eq("property_type_id", property_id).execute()
        values = [item["value"] for item in response.data] if hasattr(response, 'data') else []
        
        # Skip if no values
        if not values:
            return {
                "min": None,
                "max": None,
                "mean": None,
                "median": None,
                "std_dev": None,
                "count": 0,
                "distribution": None
            }
        
        # Filter out non-numeric values for statistics
        from statistics import mean, median, stdev
        import numpy as np
        
        numeric_values = [float(v) for v in values if isinstance(v, (int, float)) or (isinstance(v, str) and v.replace('.', '', 1).isdigit())]
        
        # Calculate statistics
        if numeric_values:
            min_val = min(numeric_values)
            max_val = max(numeric_values)
            mean_val = mean(numeric_values)
            median_val = median(numeric_values)
            std_dev_val = stdev(numeric_values) if len(numeric_values) > 1 else 0
            
            # Calculate distribution for histogram (10 bins)
            bins = 10
            hist, bin_edges = np.histogram(numeric_values, bins=bins)
            
            distribution = [
                {
                    "bin_start": float(bin_edges[i]),
                    "bin_end": float(bin_edges[i+1]),
                    "count": int(hist[i])
                }
                for i in range(bins)
            ]
        else:
            min_val = None
            max_val = None
            mean_val = None
            median_val = None
            std_dev_val = None
            distribution = None
        
        return {
            "min": min_val,
            "max": max_val,
            "mean": mean_val,
            "median": median_val,
            "std_dev": std_dev_val,
            "count": len(values),
            "distribution": distribution
        }
    except Exception as e:
        logger.error(f"Error calculating property statistics: {e}")
        return {
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "std_dev": None,
            "count": 0,
            "distribution": None
        }

if __name__ == "__main__":
    # Ensure log directory exists
    ensure_log_directory()
    
    # Run the database update
    success = perform_database_update()
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)