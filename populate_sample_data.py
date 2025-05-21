#!/usr/bin/env python3
"""
Populate the database with sample data for the CryoProtect API.
"""

import os
import sys
import logging
import psycopg2
from urllib.parse import urlparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('populate_sample_data')

# Sample molecules
SAMPLE_MOLECULES = [
    {
        "name": "Glycerol",
        "smiles": "C(C(CO)O)O",
        "pubchem_cid": "753",
        "molecular_formula": "C3H8O3",
        "molecular_weight": 92.09
    },
    {
        "name": "Dimethyl sulfoxide",
        "smiles": "CS(=O)C",
        "pubchem_cid": "679",
        "molecular_formula": "C2H6OS",
        "molecular_weight": 78.13
    },
    {
        "name": "Ethylene glycol",
        "smiles": "C(CO)O",
        "pubchem_cid": "174",
        "molecular_formula": "C2H6O2",
        "molecular_weight": 62.07
    },
    {
        "name": "Propylene glycol",
        "smiles": "CC(CO)O",
        "pubchem_cid": "1030",
        "molecular_formula": "C3H8O2",
        "molecular_weight": 76.09
    },
    {
        "name": "Sucrose",
        "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O",
        "pubchem_cid": "5988",
        "molecular_formula": "C12H22O11",
        "molecular_weight": 342.30
    },
    {
        "name": "Trehalose",
        "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O",
        "pubchem_cid": "7427",
        "molecular_formula": "C12H22O11",
        "molecular_weight": 342.30
    },
    {
        "name": "Mannitol",
        "smiles": "C(C(C(C(C(CO)O)O)O)O)O",
        "pubchem_cid": "6251",
        "molecular_formula": "C6H14O6",
        "molecular_weight": 182.17
    },
    {
        "name": "Sorbitol",
        "smiles": "C(C(C(C(C(CO)O)O)O)O)O",
        "pubchem_cid": "5780",
        "molecular_formula": "C6H14O6",
        "molecular_weight": 182.17
    },
    {
        "name": "Polyethylene glycol",
        "smiles": "C(CO)O",
        "pubchem_cid": "174",
        "molecular_formula": "C2H6O2",
        "molecular_weight": 62.07
    },
    {
        "name": "Methanol",
        "smiles": "CO",
        "pubchem_cid": "887",
        "molecular_formula": "CH4O",
        "molecular_weight": 32.04
    }
]

# Sample mixtures
SAMPLE_MIXTURES = [
    {
        "name": "Glycerol-DMSO mixture",
        "description": "Common cryoprotectant mixture"
    },
    {
        "name": "Ethylene glycol-Propylene glycol mixture",
        "description": "Automotive antifreeze mixture"
    },
    {
        "name": "Trehalose-Glycerol mixture",
        "description": "Biological cryoprotectant mixture"
    }
]

# Sample mixture components (will be filled after creating molecules and mixtures)
SAMPLE_MIXTURE_COMPONENTS = []

def get_db_connection_params():
    """Get database connection parameters from environment variables."""
    database_url = os.environ.get('DATABASE_URL')
    if not database_url:
        logger.error("DATABASE_URL environment variable not set.")
        sys.exit(1)
    
    # Parse the URL to get connection parameters
    parsed_url = urlparse(database_url)
    connection_params = {
        'host': parsed_url.hostname,
        'port': parsed_url.port or 5432,
        'dbname': parsed_url.path.lstrip('/'),
        'user': parsed_url.username,
        'password': parsed_url.password
    }
    
    logger.info(f"Database connection info: {connection_params['host']}:{connection_params['port']}/{connection_params['dbname']}")
    return connection_params

def populate_molecules(conn, cursor):
    """Populate the molecules table with sample data."""
    logger.info("Populating molecules table...")
    
    # Check if table is already populated
    cursor.execute("SELECT COUNT(*) FROM molecules")
    count = cursor.fetchone()[0]
    if count > 0:
        logger.info(f"Molecules table already has {count} records. Skipping.")
        return
    
    # Insert sample molecules
    for molecule in SAMPLE_MOLECULES:
        cursor.execute(
            """
            INSERT INTO molecules (name, smiles, pubchem_cid, molecular_formula, molecular_weight)
            VALUES (%s, %s, %s, %s, %s)
            RETURNING id
            """,
            (
                molecule["name"],
                molecule["smiles"],
                molecule["pubchem_cid"],
                molecule["molecular_formula"],
                molecule["molecular_weight"]
            )
        )
    
    conn.commit()
    logger.info(f"Inserted {len(SAMPLE_MOLECULES)} sample molecules.")

def populate_mixtures(conn, cursor):
    """Populate the mixtures table with sample data."""
    logger.info("Populating mixtures table...")
    
    # Check if table is already populated
    cursor.execute("SELECT COUNT(*) FROM mixtures")
    count = cursor.fetchone()[0]
    if count > 0:
        logger.info(f"Mixtures table already has {count} records. Skipping.")
        return
    
    # Insert sample mixtures
    for mixture in SAMPLE_MIXTURES:
        cursor.execute(
            """
            INSERT INTO mixtures (name, description)
            VALUES (%s, %s)
            RETURNING id
            """,
            (
                mixture["name"],
                mixture["description"]
            )
        )
    
    conn.commit()
    logger.info(f"Inserted {len(SAMPLE_MIXTURES)} sample mixtures.")

def populate_mixture_components(conn, cursor):
    """Populate the mixture_components table with sample data."""
    logger.info("Populating mixture_components table...")
    
    # Check if table is already populated
    cursor.execute("SELECT COUNT(*) FROM mixture_components")
    count = cursor.fetchone()[0]
    if count > 0:
        logger.info(f"Mixture_components table already has {count} records. Skipping.")
        return
    
    # Get molecule IDs
    cursor.execute("SELECT id, name FROM molecules")
    molecules = {name: id for id, name in cursor.fetchall()}
    
    # Get mixture IDs
    cursor.execute("SELECT id, name FROM mixtures")
    mixtures = {name: id for id, name in cursor.fetchall()}
    
    # Create sample mixture components
    components = [
        # Glycerol-DMSO mixture
        {"mixture_name": "Glycerol-DMSO mixture", "molecule_name": "Glycerol", "concentration": 10.0, "concentration_unit": "%v/v"},
        {"mixture_name": "Glycerol-DMSO mixture", "molecule_name": "Dimethyl sulfoxide", "concentration": 5.0, "concentration_unit": "%v/v"},
        
        # Ethylene glycol-Propylene glycol mixture
        {"mixture_name": "Ethylene glycol-Propylene glycol mixture", "molecule_name": "Ethylene glycol", "concentration": 30.0, "concentration_unit": "%v/v"},
        {"mixture_name": "Ethylene glycol-Propylene glycol mixture", "molecule_name": "Propylene glycol", "concentration": 30.0, "concentration_unit": "%v/v"},
        
        # Trehalose-Glycerol mixture
        {"mixture_name": "Trehalose-Glycerol mixture", "molecule_name": "Trehalose", "concentration": 0.5, "concentration_unit": "M"},
        {"mixture_name": "Trehalose-Glycerol mixture", "molecule_name": "Glycerol", "concentration": 10.0, "concentration_unit": "%v/v"}
    ]
    
    # Insert sample mixture components
    for component in components:
        mixture_id = mixtures.get(component["mixture_name"])
        molecule_id = molecules.get(component["molecule_name"])
        
        if not mixture_id or not molecule_id:
            logger.warning(f"Could not find mixture or molecule: {component['mixture_name']} - {component['molecule_name']}")
            continue
        
        cursor.execute(
            """
            INSERT INTO mixture_components (mixture_id, molecule_id, concentration, concentration_unit)
            VALUES (%s, %s, %s, %s)
            """,
            (
                mixture_id,
                molecule_id,
                component["concentration"],
                component["concentration_unit"]
            )
        )
    
    conn.commit()
    logger.info(f"Inserted {len(components)} sample mixture components.")

def populate_property_types(conn, cursor):
    """Populate the property_types table with sample data."""
    logger.info("Populating property_types table...")
    
    # Check if table is already populated
    cursor.execute("SELECT COUNT(*) FROM property_types")
    count = cursor.fetchone()[0]
    if count > 0:
        logger.info(f"Property_types table already has {count} records. Skipping.")
        return
    
    # Sample property types
    property_types = [
        {"name": "Melting Point", "description": "Temperature at which the solid and liquid phases coexist in equilibrium", "unit": "°C", "data_type": "numeric"},
        {"name": "Boiling Point", "description": "Temperature at which a liquid changes to a gas", "unit": "°C", "data_type": "numeric"},
        {"name": "Density", "description": "Mass per unit volume", "unit": "g/cm³", "data_type": "numeric"},
        {"name": "Solubility in Water", "description": "Maximum amount that can dissolve in water", "unit": "g/L", "data_type": "numeric"},
        {"name": "Log P", "description": "Partition coefficient between octanol and water", "unit": "", "data_type": "numeric"},
        {"name": "pKa", "description": "Acid dissociation constant", "unit": "", "data_type": "numeric"},
        {"name": "Glass Transition Temperature", "description": "Temperature at which vitrification occurs", "unit": "°C", "data_type": "numeric"},
        {"name": "Crystallization Temperature", "description": "Temperature at which crystallization occurs", "unit": "°C", "data_type": "numeric"},
        {"name": "Viscosity", "description": "Resistance to flow", "unit": "cP", "data_type": "numeric"},
        {"name": "Toxicity", "description": "Toxicity classification", "unit": "", "data_type": "text"}
    ]
    
    # Insert sample property types
    for property_type in property_types:
        cursor.execute(
            """
            INSERT INTO property_types (name, description, unit, data_type)
            VALUES (%s, %s, %s, %s)
            """,
            (
                property_type["name"],
                property_type["description"],
                property_type["unit"],
                property_type["data_type"]
            )
        )
    
    conn.commit()
    logger.info(f"Inserted {len(property_types)} sample property types.")

def populate_molecular_properties(conn, cursor):
    """Populate the molecular_properties table with sample data."""
    logger.info("Populating molecular_properties table...")
    
    # Check if table is already populated
    cursor.execute("SELECT COUNT(*) FROM molecular_properties")
    count = cursor.fetchone()[0]
    if count > 0:
        logger.info(f"Molecular_properties table already has {count} records. Skipping.")
        return
    
    # Get molecule IDs
    cursor.execute("SELECT id, name FROM molecules")
    molecules = {name: id for id, name in cursor.fetchall()}
    
    # Get property type IDs
    cursor.execute("SELECT id, name FROM property_types")
    property_types = {name: id for id, name in cursor.fetchall()}
    
    # Sample molecular properties
    properties = [
        # Glycerol
        {"molecule_name": "Glycerol", "property_name": "Melting Point", "value_numeric": 18.2, "source": "PubChem"},
        {"molecule_name": "Glycerol", "property_name": "Boiling Point", "value_numeric": 290.0, "source": "PubChem"},
        {"molecule_name": "Glycerol", "property_name": "Density", "value_numeric": 1.26, "source": "PubChem"},
        {"molecule_name": "Glycerol", "property_name": "Solubility in Water", "value_numeric": 1000.0, "source": "PubChem"},
        {"molecule_name": "Glycerol", "property_name": "Log P", "value_numeric": -1.76, "source": "PubChem"},
        
        # DMSO
        {"molecule_name": "Dimethyl sulfoxide", "property_name": "Melting Point", "value_numeric": 18.5, "source": "PubChem"},
        {"molecule_name": "Dimethyl sulfoxide", "property_name": "Boiling Point", "value_numeric": 189.0, "source": "PubChem"},
        {"molecule_name": "Dimethyl sulfoxide", "property_name": "Density", "value_numeric": 1.1, "source": "PubChem"},
        {"molecule_name": "Dimethyl sulfoxide", "property_name": "Solubility in Water", "value_numeric": 1000.0, "source": "PubChem"},
        {"molecule_name": "Dimethyl sulfoxide", "property_name": "Log P", "value_numeric": -1.35, "source": "PubChem"},
        
        # Ethylene glycol
        {"molecule_name": "Ethylene glycol", "property_name": "Melting Point", "value_numeric": -12.9, "source": "PubChem"},
        {"molecule_name": "Ethylene glycol", "property_name": "Boiling Point", "value_numeric": 197.3, "source": "PubChem"},
        {"molecule_name": "Ethylene glycol", "property_name": "Density", "value_numeric": 1.11, "source": "PubChem"},
        {"molecule_name": "Ethylene glycol", "property_name": "Solubility in Water", "value_numeric": 1000.0, "source": "PubChem"},
        {"molecule_name": "Ethylene glycol", "property_name": "Log P", "value_numeric": -1.2, "source": "PubChem"}
    ]
    
    # Insert sample molecular properties
    for prop in properties:
        molecule_id = molecules.get(prop["molecule_name"])
        property_type_id = property_types.get(prop["property_name"])
        
        if not molecule_id or not property_type_id:
            logger.warning(f"Could not find molecule or property type: {prop['molecule_name']} - {prop['property_name']}")
            continue
        
        cursor.execute(
            """
            INSERT INTO molecular_properties (molecule_id, property_type_id, value_numeric, source)
            VALUES (%s, %s, %s, %s)
            """,
            (
                molecule_id,
                property_type_id,
                prop["value_numeric"],
                prop["source"]
            )
        )
    
    conn.commit()
    logger.info(f"Inserted {len(properties)} sample molecular properties.")

def main():
    """Main function to populate the database with sample data."""
    # Get database connection parameters
    connection_params = get_db_connection_params()
    
    try:
        # Connect to the database
        conn = psycopg2.connect(**connection_params)
        cursor = conn.cursor()
        
        # Populate tables
        populate_molecules(conn, cursor)
        populate_mixtures(conn, cursor)
        populate_mixture_components(conn, cursor)
        populate_property_types(conn, cursor)
        populate_molecular_properties(conn, cursor)
        
        # Close connection
        cursor.close()
        conn.close()
        
        logger.info("Sample data population completed successfully.")
        return True
    except Exception as e:
        logger.error(f"Error populating sample data: {e}")
        return False

if __name__ == '__main__':
    success = main()
    if not success:
        sys.exit(1)
    sys.exit(0)