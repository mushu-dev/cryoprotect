#!/usr/bin/env python
"""
Convex Database Population Script

This script populates the Convex database with real scientific data,
using the migration utilities and data transformations.
"""

import os
import sys
import json
import logging
import time
from datetime import datetime
from dotenv import load_dotenv

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"convex_population_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("convex_population")

# Load environment variables
load_dotenv()

def setup_convex_client():
    """
    Setup the Convex HTTP client using the environment variables.
    """
    try:
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from database.convex_adapter import ConvexAdapter
        
        # Create client
        convex_url = os.environ.get('CONVEX_URL')
        convex_key = os.environ.get('CONVEX_DEPLOYMENT_KEY')
        
        if not convex_url:
            logger.error("CONVEX_URL environment variable not set")
            return None
            
        logger.info(f"Setting up Convex client with URL: {convex_url}")
        return ConvexAdapter(convex_url, convex_key)
    except Exception as e:
        logger.error(f"Error setting up Convex client: {str(e)}")
        return None

def fetch_molecule_data():
    """
    Fetch molecule data from the Supabase database or local files.
    """
    try:
        # Try to read from our existing data files first
        data_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                               "data", "cryoprotectant_master_list.json")
        
        if os.path.exists(data_file):
            logger.info(f"Reading molecule data from {data_file}")
            with open(data_file, 'r') as f:
                molecules = json.load(f)
                logger.info(f"Loaded {len(molecules)} molecules from file")
                return molecules
        
        # If local files don't exist, try to fetch from Supabase
        logger.info("Local data file not found, trying to fetch from Supabase")
        from supabase import create_client
        
        supabase_url = os.environ.get('SUPABASE_URL')
        supabase_key = os.environ.get('SUPABASE_KEY')
        
        if not supabase_url or not supabase_key:
            logger.error("Supabase credentials not found in environment variables")
            return []
            
        supabase = create_client(supabase_url, supabase_key)
        response = supabase.table('molecules').select('*').execute()
        
        if response.data:
            logger.info(f"Fetched {len(response.data)} molecules from Supabase")
            return response.data
        else:
            logger.warning("No molecules found in Supabase")
            return []
    
    except Exception as e:
        logger.error(f"Error fetching molecule data: {str(e)}")
        return []

def fetch_property_data(molecule_ids=None):
    """
    Fetch molecular property data from Supabase.
    
    Args:
        molecule_ids: Optional list of molecule IDs to filter by
    """
    try:
        from supabase import create_client
        
        supabase_url = os.environ.get('SUPABASE_URL')
        supabase_key = os.environ.get('SUPABASE_KEY')
        
        if not supabase_url or not supabase_key:
            logger.error("Supabase credentials not found in environment variables")
            return []
            
        supabase = create_client(supabase_url, supabase_key)
        query = supabase.table('molecular_properties').select('*')
        
        if molecule_ids:
            # Only fetch properties for specific molecules
            query = query.in_('molecule_id', molecule_ids)
            
        response = query.execute()
        
        if response.data:
            logger.info(f"Fetched {len(response.data)} molecular properties from Supabase")
            return response.data
        else:
            logger.warning("No molecular properties found in Supabase")
            return []
    
    except Exception as e:
        logger.error(f"Error fetching property data: {str(e)}")
        return []

def populate_molecules(client, molecules):
    """
    Populate the molecules table in Convex.
    
    Args:
        client: ConvexAdapter instance
        molecules: List of molecule data to populate
    
    Returns:
        dict: Mapping of original IDs to new Convex IDs
    """
    logger.info(f"Starting to populate {len(molecules)} molecules")
    id_mapping = {}
    batch_size = 50
    
    # Store start time for timing
    start_time = time.time()
    
    # Convert dictionary to list if needed
    if isinstance(molecules, dict):
        molecules_list = list(molecules.values())
    else:
        molecules_list = molecules
        
    for i in range(0, len(molecules_list), batch_size):
        batch = molecules_list[i:i+batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(len(molecules_list)-1)//batch_size + 1}, size: {len(batch)}")
        
        for molecule in batch:
            try:
                # Transform data format for Convex
                # Get the name from the molecule data
                name = "Unknown"
                if molecule.get("names") and len(molecule.get("names")) > 0:
                    name = molecule.get("names")[0]
                elif molecule.get("name"):
                    name = molecule.get("name")
                
                convex_molecule = {
                    "name": name,
                    "inchiKey": molecule.get("inchikey") or molecule.get("inchi_key", ""),
                    "smiles": molecule.get("smiles", ""),
                    "formula": molecule.get("molecular_formula") or molecule.get("formula", ""),
                    "pubchemCid": molecule.get("pubchem_cid", ""),
                    "status": molecule.get("molecule_status", "active"),
                    "consolidated": molecule.get("is_consolidated", False),
                    "createdAt": int(time.time() * 1000),
                    "updatedAt": int(time.time() * 1000)
                }
                
                # Insert into Convex
                result = client.table("molecules").insert(convex_molecule)
                
                if result.data and len(result.data) > 0:
                    # Store ID mapping
                    original_id = molecule.get("internal_id") or molecule.get("id") or molecule.get("cid") or molecule.get("pubchem_cid")
                    if original_id:
                        id_mapping[str(original_id)] = result.data[0].get("id")
                        
                    logger.debug(f"Inserted molecule {convex_molecule['name']} with ID {result.data[0].get('id')}")
                else:
                    logger.warning(f"No ID returned when inserting molecule {convex_molecule['name']}")
                
            except Exception as e:
                logger.error(f"Error inserting molecule {molecule.get('name', 'unknown')}: {str(e)}")
    
    # Calculate and log statistics
    elapsed_time = time.time() - start_time
    logger.info(f"Populated {len(id_mapping)} molecules in {elapsed_time:.2f} seconds")
    logger.info(f"Average time per molecule: {(elapsed_time / len(molecules) if molecules else 0):.4f} seconds")
    
    return id_mapping

def populate_properties(client, properties, molecule_id_mapping):
    """
    Populate the molecular_properties table in Convex.
    
    Args:
        client: ConvexAdapter instance
        properties: List of property data to populate
        molecule_id_mapping: Mapping of original molecule IDs to Convex IDs
    
    Returns:
        int: Number of properties populated
    """
    logger.info(f"Starting to populate {len(properties)} molecular properties")
    populated_count = 0
    batch_size = 100
    
    # Store start time for timing
    start_time = time.time()
    
    for i in range(0, len(properties), batch_size):
        batch = properties[i:i+batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(len(properties)-1)//batch_size + 1}, size: {len(batch)}")
        
        for prop in batch:
            try:
                # Get the Convex molecule ID
                molecule_id = molecule_id_mapping.get(str(prop.get("molecule_id")))
                
                if not molecule_id:
                    logger.warning(f"No molecule ID mapping found for property {prop.get('property_name')}, molecule_id: {prop.get('molecule_id')}")
                    continue
                
                # Transform data format for Convex
                convex_property = {
                    "moleculeId": molecule_id,
                    "propertyTypeId": str(prop.get("property_type_id", "")),
                    "propertyName": prop.get("property_name", ""),
                    "propertyType": prop.get("property_type", "custom"),
                    "numericValue": float(prop.get("numeric_value")) if prop.get("numeric_value") is not None else None,
                    "textValue": prop.get("text_value"),
                    "unit": prop.get("unit"),
                    "source": prop.get("source", "unknown"),
                    "createdAt": int(time.time() * 1000),
                    "updatedAt": int(time.time() * 1000)
                }
                
                # Insert into Convex
                result = client.table("molecularProperties").insert(convex_property)
                
                if result.data and len(result.data) > 0:
                    populated_count += 1
                    logger.debug(f"Inserted property {convex_property['propertyName']} for molecule {molecule_id}")
                else:
                    logger.warning(f"No ID returned when inserting property {convex_property['propertyName']}")
                
            except Exception as e:
                logger.error(f"Error inserting property {prop.get('property_name', 'unknown')}: {str(e)}")
    
    # Calculate and log statistics
    elapsed_time = time.time() - start_time
    logger.info(f"Populated {populated_count} properties in {elapsed_time:.2f} seconds")
    logger.info(f"Average time per property: {(elapsed_time / len(properties) if properties else 0):.4f} seconds")
    
    return populated_count

def populate_property_types(client):
    """
    Populate the property_types table in Convex with standard cryobiology properties.
    
    Args:
        client: ConvexAdapter instance
    
    Returns:
        dict: Mapping of property names to type IDs
    """
    logger.info("Populating standard property types")
    property_types = [
        {
            "name": "Glass Transition Temperature",
            "description": "Temperature at which the solution transitions to a glassy state",
            "unit": "°C",
            "dataType": "numeric",
            "category": "physical",
            "importance": "high"
        },
        {
            "name": "Melting Point",
            "description": "Temperature at which the substance melts",
            "unit": "°C",
            "dataType": "numeric",
            "category": "physical",
            "importance": "high"
        },
        {
            "name": "Molecular Weight",
            "description": "Molecular weight of the compound",
            "unit": "g/mol",
            "dataType": "numeric",
            "category": "physical",
            "importance": "high"
        },
        {
            "name": "LogP",
            "description": "Octanol-water partition coefficient",
            "unit": "",
            "dataType": "numeric",
            "category": "physical",
            "importance": "high"
        },
        {
            "name": "Solubility in Water",
            "description": "Solubility in water at room temperature",
            "unit": "g/L",
            "dataType": "numeric",
            "category": "physical",
            "importance": "high"
        },
        {
            "name": "Membrane Permeability",
            "description": "Rate at which the compound permeates cell membranes",
            "unit": "cm/s",
            "dataType": "numeric",
            "category": "biological",
            "importance": "high"
        },
        {
            "name": "Toxicity (LD50)",
            "description": "Lethal dose for 50% of test population",
            "unit": "mg/kg",
            "dataType": "numeric",
            "category": "toxicity",
            "importance": "high"
        },
        {
            "name": "Cytotoxicity (IC50)",
            "description": "Concentration causing 50% cell death",
            "unit": "μM",
            "dataType": "numeric",
            "category": "toxicity",
            "importance": "high"
        },
        {
            "name": "Osmolality at 1M",
            "description": "Osmolality of a 1M solution",
            "unit": "mOsm/kg",
            "dataType": "numeric",
            "category": "physical",
            "importance": "medium"
        },
        {
            "name": "Viscosity at 1M",
            "description": "Viscosity of a 1M solution at room temperature",
            "unit": "cP",
            "dataType": "numeric",
            "category": "physical",
            "importance": "medium"
        },
        {
            "name": "Cooling Rate Compatibility",
            "description": "Compatible cooling rates for this cryoprotectant",
            "unit": "°C/min",
            "dataType": "text",
            "category": "cryobiology",
            "importance": "high"
        },
        {
            "name": "CPA Classification",
            "description": "Classification of the cryoprotective agent",
            "unit": "",
            "dataType": "text",
            "category": "cryobiology",
            "importance": "high"
        },
        {
            "name": "Mechanism of Action",
            "description": "How the cryoprotectant works to protect cells",
            "unit": "",
            "dataType": "text",
            "category": "cryobiology",
            "importance": "medium"
        }
    ]
    
    property_name_to_id = {}
    
    for prop_type in property_types:
        try:
            # Add timestamps
            prop_type["createdAt"] = int(time.time() * 1000)
            prop_type["updatedAt"] = int(time.time() * 1000)
            
            # Insert into Convex
            result = client.table("propertyTypes").insert(prop_type)
            
            if result.data and len(result.data) > 0:
                property_name_to_id[prop_type["name"]] = result.data[0].get("id")
                logger.debug(f"Inserted property type {prop_type['name']} with ID {result.data[0].get('id')}")
            else:
                logger.warning(f"No ID returned when inserting property type {prop_type['name']}")
                
        except Exception as e:
            logger.error(f"Error inserting property type {prop_type['name']}: {str(e)}")
    
    logger.info(f"Populated {len(property_name_to_id)} property types")
    return property_name_to_id

def populate_example_experiments(client, molecule_id_mapping):
    """
    Populate the experiments table in Convex with example experiments.
    
    Args:
        client: ConvexAdapter instance
        molecule_id_mapping: Mapping of original molecule IDs to Convex IDs
    
    Returns:
        dict: Mapping of experiment names to IDs
    """
    logger.info("Populating example experiments")
    
    # First create tissue types
    tissue_types = [
        {
            "name": "Human Hepatocytes",
            "description": "Primary human liver cells",
            "species": "Human",
            "category": "Primary cells",
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        },
        {
            "name": "Mouse Embryos",
            "description": "Mouse embryos at blastocyst stage",
            "species": "Mouse",
            "category": "Embryos",
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        },
        {
            "name": "Human Oocytes",
            "description": "Human oocytes for IVF",
            "species": "Human",
            "category": "Gametes",
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        }
    ]
    
    tissue_id_mapping = {}
    for tissue in tissue_types:
        try:
            result = client.table("tissueTypes").insert(tissue)
            if result.data and len(result.data) > 0:
                tissue_id_mapping[tissue["name"]] = result.data[0].get("id")
        except Exception as e:
            logger.error(f"Error inserting tissue type {tissue['name']}: {str(e)}")
    
    # Create protocols
    protocols = [
        {
            "name": "Slow Freezing Protocol for Hepatocytes",
            "description": "Standard slow freezing protocol for human hepatocytes",
            "version": "1.0",
            "isTemplate": True,
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        },
        {
            "name": "Vitrification Protocol for Embryos",
            "description": "Vitrification protocol optimized for mouse embryos",
            "version": "2.1",
            "isTemplate": True,
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        },
        {
            "name": "Controlled-Rate Freezing for Oocytes",
            "description": "Controlled-rate freezing for human oocytes with DMSO",
            "version": "1.2",
            "isTemplate": True,
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        }
    ]
    
    protocol_id_mapping = {}
    for protocol in protocols:
        try:
            result = client.table("protocols").insert(protocol)
            if result.data and len(result.data) > 0:
                protocol_id_mapping[protocol["name"]] = result.data[0].get("id")
        except Exception as e:
            logger.error(f"Error inserting protocol {protocol['name']}: {str(e)}")
    
    # Create experiments
    current_date = datetime.now().strftime("%Y-%m-%d")
    experiments = [
        {
            "name": "Hepatocyte Slow Freezing with Glycerol",
            "description": "Testing hepatocyte viability after slow freezing with glycerol",
            "protocolId": protocol_id_mapping.get("Slow Freezing Protocol for Hepatocytes"),
            "tissueTypeId": tissue_id_mapping.get("Human Hepatocytes"),
            "experimentType": "cryopreservation",
            "startDate": current_date,
            "status": "completed",
            "researcher": "Dr. Smith",
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        },
        {
            "name": "Embryo Vitrification with DMSO/EG",
            "description": "Vitrification of mouse embryos with DMSO/EG mixture",
            "protocolId": protocol_id_mapping.get("Vitrification Protocol for Embryos"),
            "tissueTypeId": tissue_id_mapping.get("Mouse Embryos"),
            "experimentType": "cryopreservation",
            "startDate": current_date,
            "status": "completed",
            "researcher": "Dr. Jones",
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        },
        {
            "name": "Oocyte Freezing with PG",
            "description": "Controlled-rate freezing of human oocytes with propylene glycol",
            "protocolId": protocol_id_mapping.get("Controlled-Rate Freezing for Oocytes"),
            "tissueTypeId": tissue_id_mapping.get("Human Oocytes"),
            "experimentType": "cryopreservation",
            "startDate": current_date,
            "status": "in_progress",
            "researcher": "Dr. Zhang",
            "createdAt": int(time.time() * 1000),
            "updatedAt": int(time.time() * 1000)
        }
    ]
    
    experiment_id_mapping = {}
    for experiment in experiments:
        try:
            result = client.table("experiments").insert(experiment)
            if result.data and len(result.data) > 0:
                experiment_id_mapping[experiment["name"]] = result.data[0].get("id")
        except Exception as e:
            logger.error(f"Error inserting experiment {experiment['name']}: {str(e)}")
    
    # Create experiment results with real molecule IDs if available
    glycerol_id = None
    dmso_id = None
    propylene_glycol_id = None
    
    # Try to find molecule IDs by name
    for original_id, convex_id in molecule_id_mapping.items():
        molecule = client.table("molecules").select("*").eq("id", convex_id).execute()
        if molecule.data and len(molecule.data) > 0:
            mol_name = molecule.data[0].get("name", "").lower()
            if "glycerol" in mol_name:
                glycerol_id = convex_id
            elif "dimethyl sulfoxide" in mol_name or "dmso" in mol_name:
                dmso_id = convex_id
            elif "propylene glycol" in mol_name:
                propylene_glycol_id = convex_id
    
    # Create experiment results
    if experiment_id_mapping and tissue_id_mapping:
        results = [
            {
                "experimentId": experiment_id_mapping.get("Hepatocyte Slow Freezing with Glycerol"),
                "tissueTypeId": tissue_id_mapping.get("Human Hepatocytes"),
                "moleculeId": glycerol_id,
                "concentration": 10.0,
                "concentrationUnit": "%",
                "viabilityPercentage": 78.5,
                "recoveryRate": 65.2,
                "functionalityScore": 7.2,
                "timestamp": current_date,
                "createdAt": int(time.time() * 1000),
                "updatedAt": int(time.time() * 1000)
            },
            {
                "experimentId": experiment_id_mapping.get("Embryo Vitrification with DMSO/EG"),
                "tissueTypeId": tissue_id_mapping.get("Mouse Embryos"),
                "moleculeId": dmso_id,
                "concentration": 15.0,
                "concentrationUnit": "%",
                "viabilityPercentage": 92.1,
                "recoveryRate": 88.7,
                "functionalityScore": 8.9,
                "timestamp": current_date,
                "createdAt": int(time.time() * 1000),
                "updatedAt": int(time.time() * 1000)
            },
            {
                "experimentId": experiment_id_mapping.get("Oocyte Freezing with PG"),
                "tissueTypeId": tissue_id_mapping.get("Human Oocytes"),
                "moleculeId": propylene_glycol_id,
                "concentration": 12.5,
                "concentrationUnit": "%",
                "viabilityPercentage": 85.3,
                "recoveryRate": 79.8,
                "functionalityScore": 8.1,
                "timestamp": current_date,
                "createdAt": int(time.time() * 1000),
                "updatedAt": int(time.time() * 1000)
            }
        ]
        
        results_count = 0
        for result in results:
            try:
                if all([result["experimentId"], result["tissueTypeId"]]):
                    res = client.table("experimentResults").insert(result)
                    if res.data and len(res.data) > 0:
                        results_count += 1
            except Exception as e:
                logger.error(f"Error inserting experiment result: {str(e)}")
        
        logger.info(f"Populated {results_count} experiment results")
    
    return experiment_id_mapping

def main():
    """
    Main function to run the Convex database population.
    """
    logger.info("Starting Convex database population")
    
    # Setup Convex client
    client = setup_convex_client()
    if not client:
        logger.error("Failed to setup Convex client, exiting")
        return
    
    # Populate property types first
    property_types = populate_property_types(client)
    logger.info(f"Populated {len(property_types)} property types")
    
    # Fetch and populate molecules
    molecules = fetch_molecule_data()
    if not molecules:
        logger.error("No molecule data found, exiting")
        return
    
    # Limit to a smaller set for initial testing
    test_limit = int(os.environ.get('TEST_LIMIT', '0'))
    if test_limit > 0:
        logger.info(f"Using test limit of {test_limit} molecules")
        molecules = molecules[:test_limit]
    
    # Populate molecules
    molecule_id_mapping = populate_molecules(client, molecules)
    logger.info(f"Populated {len(molecule_id_mapping)} molecules with ID mapping")
    
    # Fetch and populate properties
    properties = fetch_property_data(list(molecule_id_mapping.keys()))
    if properties:
        property_count = populate_properties(client, properties, molecule_id_mapping)
        logger.info(f"Populated {property_count} molecular properties")
    else:
        logger.warning("No property data found, skipping property population")
    
    # Populate example experiments
    experiment_id_mapping = populate_example_experiments(client, molecule_id_mapping)
    logger.info(f"Populated {len(experiment_id_mapping)} example experiments")
    
    logger.info("Convex database population completed")

if __name__ == "__main__":
    main()