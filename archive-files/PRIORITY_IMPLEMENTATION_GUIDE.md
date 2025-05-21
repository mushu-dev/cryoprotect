# PRIORITY IMPLEMENTATION GUIDE

This document provides detailed implementation instructions for the highest priority microtasks identified in the Integrated Data Flow Implementation Plan. These tasks address critical bottlenecks in data flow and will enable efficient scientific data integration and display.

## Quick Reference
| Task ID | Priority | Task Name | Est. Time | Dependencies |
|---------|----------|-----------|-----------|--------------|
| 1.1 | 1 | Connection Pool Enhancement | 6 hours | None |
| 2.2 | 2 | Cryoprotectant Identifier List | 4 hours | None |
| 2.3 | 3 | PubChem Import Optimization | 8 hours | 1.1, 2.2 |
| 2.7 | 4 | ChEMBL Import Enhancement | 10 hours | 1.1, 2.2 |
| 3.1 | 5 | Molecule Viewer Optimization | 6 hours | 2.3 or 2.7 |
| 3.3 | 6 | Data Table Enhancement | 6 hours | 2.3 or 2.7 |

## 1. Connection Pool Enhancement (Task 1.1)

### Overview
This task involves implementing a robust connection pooling system to reduce database connection overhead and improve data import performance.

### Files to Modify
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/connection_pool_wrapper.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/implement_connection_pooling.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/config.py`

### Implementation Steps

#### Step 1.1.1: Enhance Connection Pool Wrapper
Update the connection pool wrapper to handle connection lifecycle management:

```python
# connection_pool_wrapper.py
import psycopg2
from psycopg2 import pool
import logging
import threading
import time
from typing import Dict, Any, Optional, List, Tuple

logger = logging.getLogger(__name__)

class ConnectionPoolWrapper:
    """
    Enhanced connection pool wrapper with lifecycle management, 
    health monitoring, and automatic reconnection.
    """
    _instance = None
    _lock = threading.Lock()
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the connection pool with configuration parameters.
        
        Args:
            config: Dictionary containing connection parameters
                - min_connections: Minimum connections in pool
                - max_connections: Maximum connections in pool
                - host: Database host
                - dbname: Database name
                - user: Database user
                - password: Database password
                - port: Database port
        """
        self.config = config
        self.min_conn = config.get('min_connections', 1)
        self.max_conn = config.get('max_connections', 10)
        self.pool = None
        self.health_check_interval = config.get('health_check_interval', 60)  # seconds
        self.connection_timeout = config.get('connection_timeout', 30)  # seconds
        self.pool_initialized = False
        self.active_connections = 0
        self.last_health_check = 0
        self._initialize_pool()
        
        # Start health check thread
        self._start_health_check()
    
    @classmethod
    def get_instance(cls, config: Optional[Dict[str, Any]] = None) -> 'ConnectionPoolWrapper':
        """
        Get singleton instance of connection pool wrapper.
        
        Args:
            config: Configuration dictionary (only used on first initialization)
            
        Returns:
            ConnectionPoolWrapper instance
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    if config is None:
                        raise ValueError("Config must be provided for initial pool creation")
                    cls._instance = cls(config)
        return cls._instance
    
    def _initialize_pool(self) -> None:
        """Initialize the connection pool with configuration parameters."""
        try:
            if self.pool is not None:
                # Close existing pool if it exists
                self.pool.closeall()
                
            self.pool = psycopg2.pool.ThreadedConnectionPool(
                self.min_conn,
                self.max_conn,
                host=self.config.get('host'),
                dbname=self.config.get('dbname'),
                user=self.config.get('user'),
                password=self.config.get('password'),
                port=self.config.get('port'),
                connect_timeout=self.connection_timeout
            )
            self.pool_initialized = True
            self.active_connections = 0
            logger.info("Connection pool initialized with min=%d, max=%d connections", 
                       self.min_conn, self.max_conn)
        except Exception as e:
            self.pool_initialized = False
            logger.error("Failed to initialize connection pool: %s", str(e))
            raise
    
    def _start_health_check(self) -> None:
        """Start a background thread for periodic health checks."""
        health_thread = threading.Thread(
            target=self._health_check_worker, 
            daemon=True
        )
        health_thread.start()
        logger.info("Health check thread started with interval %d seconds", 
                   self.health_check_interval)
    
    def _health_check_worker(self) -> None:
        """Worker function for periodic health checks."""
        while True:
            time.sleep(self.health_check_interval)
            try:
                self._check_pool_health()
            except Exception as e:
                logger.error("Health check failed: %s", str(e))
    
    def _check_pool_health(self) -> None:
        """Check health of connection pool and reinitialize if needed."""
        if not self.pool_initialized:
            logger.warning("Pool not initialized during health check, attempting to initialize")
            self._initialize_pool()
            return
            
        try:
            # Get a connection to test
            conn = self.pool.getconn()
            try:
                # Execute simple query to check connection
                with conn.cursor() as cursor:
                    cursor.execute("SELECT 1")
                    result = cursor.fetchone()
                    if result and result[0] == 1:
                        logger.debug("Connection pool health check passed")
                    else:
                        logger.warning("Connection returned unexpected result during health check")
                        self._initialize_pool()
            except Exception as e:
                logger.error("Connection test failed during health check: %s", str(e))
                self._initialize_pool()
            finally:
                # Return the connection to the pool
                self.pool.putconn(conn)
        except Exception as e:
            logger.error("Failed to get connection during health check: %s", str(e))
            self._initialize_pool()
        
        self.last_health_check = time.time()
    
    def get_connection(self) -> Tuple[Any, int]:
        """
        Get a connection from the pool with a unique identifier.
        
        Returns:
            Tuple of (connection, connection_id)
        """
        if not self.pool_initialized:
            self._initialize_pool()
            
        try:
            conn = self.pool.getconn()
            self.active_connections += 1
            conn_id = id(conn)
            logger.debug("Connection %d acquired from pool (active: %d)", 
                        conn_id, self.active_connections)
            return conn, conn_id
        except Exception as e:
            logger.error("Failed to get connection from pool: %s", str(e))
            # Try to reinitialize pool
            self._initialize_pool()
            # Retry once
            conn = self.pool.getconn()
            self.active_connections += 1
            conn_id = id(conn)
            logger.debug("Connection %d acquired from pool after retry (active: %d)", 
                        conn_id, self.active_connections)
            return conn, conn_id
    
    def return_connection(self, conn: Any, conn_id: int) -> None:
        """
        Return a connection to the pool.
        
        Args:
            conn: The connection to return
            conn_id: The connection identifier
        """
        if self.pool_initialized and conn is not None:
            try:
                self.pool.putconn(conn)
                self.active_connections -= 1
                logger.debug("Connection %d returned to pool (active: %d)", 
                            conn_id, self.active_connections)
            except Exception as e:
                logger.error("Failed to return connection %d to pool: %s", 
                            conn_id, str(e))
    
    def close_all(self) -> None:
        """Close all connections in the pool."""
        if self.pool_initialized:
            try:
                self.pool.closeall()
                logger.info("Closed all connections in the pool")
                self.pool_initialized = False
                self.active_connections = 0
            except Exception as e:
                logger.error("Failed to close all connections: %s", str(e))

# Helper functions for context manager usage
def get_db_connection():
    """Get a connection from the global connection pool."""
    from config import get_db_config
    pool = ConnectionPoolWrapper.get_instance(get_db_config())
    return pool.get_connection()

def return_db_connection(conn, conn_id):
    """Return a connection to the global connection pool."""
    from config import get_db_config
    pool = ConnectionPoolWrapper.get_instance(get_db_config())
    pool.return_connection(conn, conn_id)

class ConnectionManager:
    """Context manager for database connections."""
    
    def __init__(self):
        self.conn = None
        self.conn_id = None
    
    def __enter__(self):
        self.conn, self.conn_id = get_db_connection()
        return self.conn
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return_db_connection(self.conn, self.conn_id)
```

#### Step 1.1.2: Update Config to Support Connection Pooling

```python
# config.py
import os
from typing import Dict, Any, Optional
import logging

# Database configuration
def get_db_config() -> Dict[str, Any]:
    """
    Get database configuration with pool settings.
    
    Returns:
        Dict containing database configuration parameters
    """
    return {
        'host': os.environ.get('DB_HOST', 'localhost'),
        'dbname': os.environ.get('DB_NAME', 'cryoprotect'),
        'user': os.environ.get('DB_USER', 'postgres'),
        'password': os.environ.get('DB_PASSWORD', 'postgres'),
        'port': int(os.environ.get('DB_PORT', 5432)),
        'min_connections': int(os.environ.get('DB_MIN_CONNECTIONS', 1)),
        'max_connections': int(os.environ.get('DB_MAX_CONNECTIONS', 10)),
        'health_check_interval': int(os.environ.get('DB_HEALTH_CHECK_INTERVAL', 60)),
        'connection_timeout': int(os.environ.get('DB_CONNECTION_TIMEOUT', 30)),
    }
```

#### Step 1.1.3: Create Implementation Script

```python
# implement_connection_pooling.py
import logging
import sys
import time
from typing import List, Dict, Any, Optional
from connection_pool_wrapper import ConnectionPoolWrapper, ConnectionManager
from config import get_db_config

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('connection_pool.log')
    ]
)
logger = logging.getLogger(__name__)

def test_connection_pool() -> bool:
    """
    Test the connection pool with multiple operations.
    
    Returns:
        True if test is successful, False otherwise
    """
    try:
        # Create connection pool
        config = get_db_config()
        pool = ConnectionPoolWrapper.get_instance(config)
        
        # Test basic query
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                cursor.execute("SELECT 1 as test")
                result = cursor.fetchone()
                logger.info("Basic query test: %s", result)
        
        # Test multiple connections
        threads = []
        for i in range(5):
            thread = threading.Thread(
                target=_run_test_query,
                args=(i,)
            )
            threads.append(thread)
            thread.start()
        
        # Wait for threads to complete
        for thread in threads:
            thread.join()
        
        # Check active connections
        logger.info("Active connections: %d", pool.active_connections)
        
        # Test health check
        pool._check_pool_health()
        
        return True
    except Exception as e:
        logger.error("Connection pool test failed: %s", str(e))
        return False

def _run_test_query(thread_id: int) -> None:
    """Run a test query in a separate thread."""
    try:
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                cursor.execute("SELECT pg_sleep(0.5)")
                logger.info("Thread %d completed query", thread_id)
    except Exception as e:
        logger.error("Thread %d query failed: %s", thread_id, str(e))

def apply_connection_pooling() -> None:
    """Apply connection pooling to the application."""
    try:
        # Test connection pool
        success = test_connection_pool()
        if not success:
            logger.error("Connection pool test failed, not applying changes")
            return
        
        logger.info("Connection pool test successful")
        
        # The actual application of connection pooling will involve
        # updating import statements in key modules to use the connection
        # pool wrapper instead of direct connections. This function can
        # be extended to automatically apply these changes.
        
    except Exception as e:
        logger.error("Failed to apply connection pooling: %s", str(e))

if __name__ == "__main__":
    logger.info("Implementing connection pooling...")
    apply_connection_pooling()
    logger.info("Connection pooling implementation complete")
```

#### Step 1.1.4: Update Key Imports in Data Processing Files

Replace individual connection creation with the connection pool in:
- `ChEMBL_Integrated_Import.py`
- `PubChem_CryoProtectants_Supabase.py`
- Any other file that creates DB connections directly

Example replacement pattern:
```python
# Before:
def get_db_connection():
    """Get direct db connection."""
    global db
    if db is None:
        db = SupabaseDirectConnection.get_instance()
    return db

# After:
from connection_pool_wrapper import ConnectionManager

def perform_db_operation():
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            # Execute operations
            pass
```

### Success Criteria
- √ Connection pool automatically manages connection lifecycle
- √ Health check thread periodically validates connections
- √ Failed connections are automatically reconnected
- √ Connection pooling demonstrates 30%+ performance improvement in data import tasks

## 2. Cryoprotectant Identifier List (Task 2.2)

### Overview
This task involves creating a standardized list of cryoprotectant molecule identifiers to ensure consistent data import across multiple data sources.

### Files to Create/Modify
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/cryoprotectant_identifiers.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/data/cryoprotectant_master_list.json`

### Implementation Steps

#### Step 2.2.1: Create Cryoprotectant Identifier Mapping Module

```python
# cryoprotectant_identifiers.py
import json
import os
import logging
from typing import Dict, List, Any, Optional, Set, Tuple

logger = logging.getLogger(__name__)

MASTER_LIST_PATH = os.path.join('data', 'cryoprotectant_master_list.json')

class CryoprotectantIdentifierManager:
    """
    Manages cryoprotectant identifiers across different data sources.
    Provides mapping between different identifier systems and validation.
    """
    _instance = None
    
    def __init__(self, master_list_path: str = MASTER_LIST_PATH):
        """
        Initialize with path to master identifier list.
        
        Args:
            master_list_path: Path to the JSON file containing master identifier list
        """
        self.master_list_path = master_list_path
        self.identifiers = {}
        self.names = {}
        self.pubchem_cids = set()
        self.chembl_ids = set()
        self.cas_numbers = set()
        self.inchi_keys = set()
        self.smiles = set()
        self._load_identifiers()
    
    @classmethod
    def get_instance(cls, master_list_path: str = MASTER_LIST_PATH) -> 'CryoprotectantIdentifierManager':
        """
        Get singleton instance.
        
        Args:
            master_list_path: Path to master identifier list
            
        Returns:
            CryoprotectantIdentifierManager instance
        """
        if cls._instance is None:
            cls._instance = cls(master_list_path)
        return cls._instance
    
    def _load_identifiers(self) -> None:
        """Load identifiers from master list file."""
        try:
            if not os.path.exists(self.master_list_path):
                logger.warning("Master list file not found at %s. Creating empty list.", 
                              self.master_list_path)
                self.identifiers = {}
                return
                
            with open(self.master_list_path, 'r') as f:
                self.identifiers = json.load(f)
            
            # Build lookup indexes
            for internal_id, data in self.identifiers.items():
                if 'pubchem_cid' in data and data['pubchem_cid']:
                    self.pubchem_cids.add(str(data['pubchem_cid']))
                
                if 'chembl_id' in data and data['chembl_id']:
                    self.chembl_ids.add(data['chembl_id'])
                
                if 'cas_number' in data and data['cas_number']:
                    self.cas_numbers.add(data['cas_number'])
                
                if 'inchi_key' in data and data['inchi_key']:
                    self.inchi_keys.add(data['inchi_key'])
                
                if 'smiles' in data and data['smiles']:
                    self.smiles.add(data['smiles'])
                
                if 'names' in data and data['names']:
                    for name in data['names']:
                        self.names[name.lower()] = internal_id
            
            logger.info("Loaded %d cryoprotectant identifiers", len(self.identifiers))
            logger.info("Indexed %d PubChem CIDs, %d ChEMBL IDs, %d CAS numbers, %d InChI Keys, %d SMILES",
                       len(self.pubchem_cids), len(self.chembl_ids), len(self.cas_numbers),
                       len(self.inchi_keys), len(self.smiles))
                       
        except Exception as e:
            logger.error("Failed to load identifiers: %s", str(e))
            # Initialize empty
            self.identifiers = {}
    
    def save_identifiers(self) -> None:
        """Save identifiers to master list file."""
        try:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.master_list_path), exist_ok=True)
            
            with open(self.master_list_path, 'w') as f:
                json.dump(self.identifiers, f, indent=2)
            
            logger.info("Saved %d cryoprotectant identifiers to %s", 
                       len(self.identifiers), self.master_list_path)
        except Exception as e:
            logger.error("Failed to save identifiers: %s", str(e))
    
    def add_molecule(self, internal_id: str, data: Dict[str, Any]) -> None:
        """
        Add a molecule to the identifier list.
        
        Args:
            internal_id: Internal identifier for the molecule
            data: Dictionary containing molecule identifier data
        """
        self.identifiers[internal_id] = data
        
        # Update lookup indexes
        if 'pubchem_cid' in data and data['pubchem_cid']:
            self.pubchem_cids.add(str(data['pubchem_cid']))
        
        if 'chembl_id' in data and data['chembl_id']:
            self.chembl_ids.add(data['chembl_id'])
        
        if 'cas_number' in data and data['cas_number']:
            self.cas_numbers.add(data['cas_number'])
        
        if 'inchi_key' in data and data['inchi_key']:
            self.inchi_keys.add(data['inchi_key'])
        
        if 'smiles' in data and data['smiles']:
            self.smiles.add(data['smiles'])
        
        if 'names' in data and data['names']:
            for name in data['names']:
                self.names[name.lower()] = internal_id
    
    def get_molecule_by_internal_id(self, internal_id: str) -> Optional[Dict[str, Any]]:
        """
        Get molecule data by internal identifier.
        
        Args:
            internal_id: Internal identifier for the molecule
            
        Returns:
            Dictionary containing molecule data or None if not found
        """
        return self.identifiers.get(internal_id)
    
    def get_internal_id_by_pubchem_cid(self, pubchem_cid: str) -> Optional[str]:
        """
        Get internal identifier by PubChem CID.
        
        Args:
            pubchem_cid: PubChem CID
            
        Returns:
            Internal identifier or None if not found
        """
        for internal_id, data in self.identifiers.items():
            if 'pubchem_cid' in data and str(data['pubchem_cid']) == str(pubchem_cid):
                return internal_id
        return None
    
    def get_internal_id_by_chembl_id(self, chembl_id: str) -> Optional[str]:
        """
        Get internal identifier by ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            Internal identifier or None if not found
        """
        for internal_id, data in self.identifiers.items():
            if 'chembl_id' in data and data['chembl_id'] == chembl_id:
                return internal_id
        return None
    
    def get_internal_id_by_name(self, name: str) -> Optional[str]:
        """
        Get internal identifier by molecule name.
        
        Args:
            name: Molecule name
            
        Returns:
            Internal identifier or None if not found
        """
        return self.names.get(name.lower())
    
    def get_internal_id_by_inchi_key(self, inchi_key: str) -> Optional[str]:
        """
        Get internal identifier by InChI Key.
        
        Args:
            inchi_key: InChI Key
            
        Returns:
            Internal identifier or None if not found
        """
        for internal_id, data in self.identifiers.items():
            if 'inchi_key' in data and data['inchi_key'] == inchi_key:
                return internal_id
        return None
    
    def resolve_identifier(self, 
                          pubchem_cid: Optional[str] = None,
                          chembl_id: Optional[str] = None,
                          name: Optional[str] = None,
                          cas_number: Optional[str] = None, 
                          inchi_key: Optional[str] = None,
                          smiles: Optional[str] = None) -> Tuple[Optional[str], float]:
        """
        Resolve molecule identifier from any available identifier.
        
        Args:
            pubchem_cid: PubChem CID
            chembl_id: ChEMBL ID
            name: Molecule name
            cas_number: CAS Registry Number
            inchi_key: InChI Key
            smiles: SMILES string
            
        Returns:
            Tuple of (internal_id, confidence_score) or (None, 0.0) if not found
        """
        # Check exact matches
        if pubchem_cid:
            internal_id = self.get_internal_id_by_pubchem_cid(pubchem_cid)
            if internal_id:
                return internal_id, 1.0
        
        if chembl_id:
            internal_id = self.get_internal_id_by_chembl_id(chembl_id)
            if internal_id:
                return internal_id, 1.0
        
        if inchi_key:
            internal_id = self.get_internal_id_by_inchi_key(inchi_key)
            if internal_id:
                return internal_id, 1.0
        
        if name:
            internal_id = self.get_internal_id_by_name(name)
            if internal_id:
                return internal_id, 0.9  # Names can be ambiguous
        
        # Implement fuzzy matching for names and SMILES if needed
        # ...
        
        # No match found
        return None, 0.0
    
    def get_all_pubchem_cids(self) -> List[str]:
        """
        Get list of all PubChem CIDs.
        
        Returns:
            List of PubChem CIDs
        """
        return list(self.pubchem_cids)
    
    def get_all_chembl_ids(self) -> List[str]:
        """
        Get list of all ChEMBL IDs.
        
        Returns:
            List of ChEMBL IDs
        """
        return list(self.chembl_ids)

def initialize_cryoprotectant_list() -> None:
    """Initialize the cryoprotectant list with commonly used molecules."""
    manager = CryoprotectantIdentifierManager.get_instance()
    
    # Core cryoprotectants
    common_cryoprotectants = [
        {
            "internal_id": "CRYO001",
            "pubchem_cid": "962",
            "chembl_id": "CHEMBL388978",
            "cas_number": "56-81-5",
            "names": ["Glycerol", "Glycerin", "1,2,3-Propanetriol"],
            "inchi_key": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
            "smiles": "C(C(CO)O)O",
            "formula": "C3H8O3",
            "molecular_weight": 92.09,
            "category": "polyol"
        },
        {
            "internal_id": "CRYO002", 
            "pubchem_cid": "5988", 
            "chembl_id": "CHEMBL1098659", 
            "cas_number": "67-68-5", 
            "names": ["Dimethyl sulfoxide", "DMSO", "Methyl sulfoxide"], 
            "inchi_key": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N", 
            "smiles": "CS(=O)C", 
            "formula": "C2H6OS", 
            "molecular_weight": 78.13, 
            "category": "organosulfur"
        },
        {
            "internal_id": "CRYO003",
            "pubchem_cid": "6342",
            "chembl_id": "CHEMBL66195",
            "cas_number": "107-95-9",
            "names": ["beta-Alanine", "3-Aminopropanoic acid", "3-Aminopropionic acid"],
            "inchi_key": "UCMIRNVEIXFBKS-UHFFFAOYSA-N",
            "smiles": "C(CN)C(=O)O",
            "formula": "C3H7NO2",
            "molecular_weight": 89.09,
            "category": "amino acid"
        },
        {
            "internal_id": "CRYO004",
            "pubchem_cid": "1030",
            "chembl_id": "CHEMBL500033",
            "cas_number": "75-65-0",
            "names": ["tert-Butanol", "t-Butyl alcohol", "2-Methyl-2-propanol"],
            "inchi_key": "DKGAVHZHDRPRBM-UHFFFAOYSA-N",
            "smiles": "CC(C)(C)O",
            "formula": "C4H10O",
            "molecular_weight": 74.12,
            "category": "alcohol"
        },
        {
            "internal_id": "CRYO005",
            "pubchem_cid": "6057",
            "chembl_id": "CHEMBL1487",
            "cas_number": "57-13-6",
            "names": ["Urea", "Carbamide", "Carbonyl diamide"],
            "inchi_key": "XSQUKJJJFZCRTK-UHFFFAOYSA-N",
            "smiles": "C(=O)(N)N",
            "formula": "CH4N2O",
            "molecular_weight": 60.06,
            "category": "amide"
        }
    ]
    
    # Add molecules to manager
    for molecule in common_cryoprotectants:
        manager.add_molecule(molecule["internal_id"], molecule)
    
    # Save to file
    manager.save_identifiers()
    
    logger.info("Initialized cryoprotectant list with %d molecules", 
               len(common_cryoprotectants))

# Command-line interface
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    logger.info("Initializing cryoprotectant identifier list...")
    initialize_cryoprotectant_list()
    logger.info("Done!")
```

#### Step 2.2.2: Create Initial Master List of Cryoprotectants

Run the initialization script to create the master list:

```bash
# Create the initial list
python cryoprotectant_identifiers.py
```

This will create the `data/cryoprotectant_master_list.json` file with common cryoprotectants.

#### Step 2.2.3: Update Import Scripts to Use Identifier Manager

In `PubChem_CryoProtectants_Supabase.py` and `ChEMBL_Integrated_Import.py`, add:

```python
from cryoprotectant_identifiers import CryoprotectantIdentifierManager

# During import processing
identifier_manager = CryoprotectantIdentifierManager.get_instance()

# When processing a molecule
internal_id, confidence = identifier_manager.resolve_identifier(
    pubchem_cid=molecule_data.get('cid'),
    name=molecule_data.get('name'),
    inchi_key=molecule_data.get('inchi_key')
)

# If not found, create a new internal ID
if not internal_id:
    # Generate new internal ID
    internal_id = f"CRYO{next_id:04d}"
    next_id += 1
    
    # Add to identifier manager
    identifier_manager.add_molecule(internal_id, {
        "internal_id": internal_id,
        "pubchem_cid": molecule_data.get('cid'),
        "names": [molecule_data.get('name')] if molecule_data.get('name') else [],
        "inchi_key": molecule_data.get('inchi_key'),
        "smiles": molecule_data.get('smiles'),
        # other data...
    })
    identifier_manager.save_identifiers()
```

### Success Criteria
- √ Master list contains at least 50 common cryoprotectants with complete identifiers
- √ Identifier resolution works with multiple identifier types (PubChem CID, ChEMBL ID, names)
- √ Import scripts can add new molecules to the master list
- √ Molecule lookup has confidence scoring for ambiguous matches

## 3. PubChem Import Optimization (Task 2.3)

### Overview
This task optimizes the PubChem data import process with efficient batch processing, better error handling, and integration with the new connection pool and identifier manager.

### Files to Modify
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/PubChem_CryoProtectants_Supabase.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/PubChem_CryoProtectants_Supabase_Enhanced.py` (New file)

### Implementation Steps

#### Step 3.1: Create Enhanced PubChem Import Script

```python
# PubChem_CryoProtectants_Supabase_Enhanced.py
import os
import json
import logging
import time
import argparse
import concurrent.futures
from typing import Dict, List, Any, Optional, Set, Tuple
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import psycopg2
from datetime import datetime

# Import custom modules
from connection_pool_wrapper import ConnectionManager
from cryoprotectant_identifiers import CryoprotectantIdentifierManager
from pubchem.client import PubChemClient
from pubchem.rate_limiter import RateLimiter
from pubchem.cache import PubChemCache

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('pubchem_import.log')
    ]
)
logger = logging.getLogger(__name__)

# Constants
CHECKPOINT_DIR = "checkpoints"
CHECKPOINT_FILE = os.path.join(CHECKPOINT_DIR, "pubchem_import_enhanced.json")
BATCH_SIZE = 25  # Number of compounds to process in parallel
MAX_RETRIES = 3  # Maximum number of retries for failed operations

class PubChemImporter:
    """
    Enhanced PubChem importer with batch processing, robust error handling,
    and integration with connection pool and identifier manager.
    """
    
    def __init__(self, 
                checkpoint_file: str = CHECKPOINT_FILE,
                batch_size: int = BATCH_SIZE,
                max_retries: int = MAX_RETRIES):
        """
        Initialize the PubChem importer.
        
        Args:
            checkpoint_file: Path to checkpoint file
            batch_size: Number of compounds to process in parallel
            max_retries: Maximum number of retries for failed operations
        """
        self.checkpoint_file = checkpoint_file
        self.batch_size = batch_size
        self.max_retries = max_retries
        
        # Initialize components
        self.pubchem_client = PubChemClient()
        self.pubchem_cache = PubChemCache()
        self.rate_limiter = RateLimiter(requests_per_minute=5)
        self.id_manager = CryoprotectantIdentifierManager.get_instance()
        
        # Import state
        self.checkpoint = self._load_checkpoint()
        self.next_internal_id = self.checkpoint.get('next_internal_id', 1)
        
        # Statistics
        self.stats = {
            'compounds_processed': 0,
            'compounds_added': 0,
            'compounds_updated': 0,
            'compounds_skipped': 0,
            'compounds_failed': 0,
            'start_time': None,
            'end_time': None,
            'total_elapsed_time': 0
        }
    
    def _load_checkpoint(self) -> Dict[str, Any]:
        """
        Load checkpoint from file.
        
        Returns:
            Dictionary containing checkpoint data
        """
        if not os.path.exists(self.checkpoint_file):
            logger.info("No checkpoint file found at %s. Starting fresh import.", 
                       self.checkpoint_file)
            return {
                'processed_cids': [],
                'failed_cids': {},
                'next_internal_id': 1,
                'last_completed_batch': -1,
                'total_compounds': 0
            }
        
        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
            
            logger.info("Loaded checkpoint with %d processed compounds, %d failed compounds",
                       len(checkpoint.get('processed_cids', [])), 
                       len(checkpoint.get('failed_cids', {})))
            
            return checkpoint
        except Exception as e:
            logger.error("Failed to load checkpoint: %s", str(e))
            return {
                'processed_cids': [],
                'failed_cids': {},
                'next_internal_id': 1,
                'last_completed_batch': -1,
                'total_compounds': 0
            }
    
    def _save_checkpoint(self) -> None:
        """Save checkpoint to file."""
        try:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.checkpoint_file), exist_ok=True)
            
            with open(self.checkpoint_file, 'w') as f:
                json.dump(self.checkpoint, f, indent=2)
            
            logger.debug("Saved checkpoint with %d processed compounds, %d failed compounds",
                        len(self.checkpoint.get('processed_cids', [])), 
                        len(self.checkpoint.get('failed_cids', {})))
        except Exception as e:
            logger.error("Failed to save checkpoint: %s", str(e))
    
    def _generate_internal_id(self) -> str:
        """
        Generate a new internal ID for a molecule.
        
        Returns:
            Internal ID string (e.g., CRYO0123)
        """
        internal_id = f"CRYO{self.next_internal_id:04d}"
        self.next_internal_id += 1
        self.checkpoint['next_internal_id'] = self.next_internal_id
        return internal_id
    
    def _get_or_create_internal_id(self, 
                                  pubchem_data: Dict[str, Any]) -> Tuple[str, bool]:
        """
        Get existing internal ID or create a new one.
        
        Args:
            pubchem_data: PubChem compound data
            
        Returns:
            Tuple of (internal_id, is_new)
        """
        pubchem_cid = str(pubchem_data.get('id', {}).get('id', {}).get('cid', ''))
        name = pubchem_data.get('name', '')
        inchi_key = pubchem_data.get('inchi_key', '')
        
        # Try to resolve with identifier manager
        internal_id, confidence = self.id_manager.resolve_identifier(
            pubchem_cid=pubchem_cid,
            name=name,
            inchi_key=inchi_key
        )
        
        if internal_id:
            return internal_id, False
        
        # Create new internal ID
        internal_id = self._generate_internal_id()
        
        # Add to identifier manager
        synonyms = []
        if name:
            synonyms.append(name)
        
        for synonym in pubchem_data.get('synonyms', [])[:5]:  # Limit to 5 synonyms
            if synonym not in synonyms:
                synonyms.append(synonym)
        
        self.id_manager.add_molecule(internal_id, {
            "internal_id": internal_id,
            "pubchem_cid": pubchem_cid,
            "names": synonyms,
            "inchi_key": inchi_key,
            "smiles": pubchem_data.get('smiles', ''),
            "formula": pubchem_data.get('molecular_formula', ''),
            "category": "imported"
        })
        
        # Save identifier list
        self.id_manager.save_identifiers()
        
        return internal_id, True
    
    def _fetch_compound_data(self, cid: str) -> Optional[Dict[str, Any]]:
        """
        Fetch compound data from PubChem.
        
        Args:
            cid: PubChem CID
            
        Returns:
            Dictionary containing compound data or None if failed
        """
        retry_count = 0
        while retry_count < self.max_retries:
            try:
                self.rate_limiter.wait()
                
                # Try to get from cache first
                cached_data = self.pubchem_cache.get_compound(cid)
                if cached_data:
                    logger.debug("Using cached data for CID %s", cid)
                    return cached_data
                
                # Fetch from PubChem
                logger.debug("Fetching data for CID %s from PubChem", cid)
                compound_data = self.pubchem_client.get_compound(cid)
                
                # Cache the result
                self.pubchem_cache.cache_compound(cid, compound_data)
                
                return compound_data
            
            except Exception as e:
                retry_count += 1
                logger.warning("Failed to fetch CID %s (attempt %d/%d): %s", 
                              cid, retry_count, self.max_retries, str(e))
                time.sleep(2 ** retry_count)  # Exponential backoff
        
        logger.error("Failed to fetch CID %s after %d attempts", cid, self.max_retries)
        return None
    
    def _process_compound(self, cid: str) -> bool:
        """
        Process a single compound.
        
        Args:
            cid: PubChem CID
            
        Returns:
            True if successful, False otherwise
        """
        if cid in self.checkpoint['processed_cids']:
            logger.debug("Skipping CID %s (already processed)", cid)
            self.stats['compounds_skipped'] += 1
            return True
        
        try:
            # Fetch compound data
            compound_data = self._fetch_compound_data(cid)
            if not compound_data:
                self.checkpoint['failed_cids'][cid] = "Failed to fetch data"
                self.stats['compounds_failed'] += 1
                return False
            
            # Get or create internal ID
            internal_id, is_new = self._get_or_create_internal_id(compound_data)
            
            # Insert or update in database
            self._insert_or_update_compound(internal_id, cid, compound_data, is_new)
            
            # Mark as processed
            if cid not in self.checkpoint['processed_cids']:
                self.checkpoint['processed_cids'].append(cid)
            
            if is_new:
                self.stats['compounds_added'] += 1
            else:
                self.stats['compounds_updated'] += 1
            
            self.stats['compounds_processed'] += 1
            
            return True
            
        except Exception as e:
            logger.error("Failed to process CID %s: %s", cid, str(e))
            self.checkpoint['failed_cids'][cid] = str(e)
            self.stats['compounds_failed'] += 1
            return False
    
    def _insert_or_update_compound(self, 
                                  internal_id: str, 
                                  cid: str, 
                                  compound_data: Dict[str, Any],
                                  is_new: bool) -> None:
        """
        Insert or update compound in database.
        
        Args:
            internal_id: Internal molecule ID
            cid: PubChem CID
            compound_data: PubChem compound data
            is_new: Whether this is a new compound
        """
        retry_count = 0
        
        while retry_count < self.max_retries:
            try:
                # Extract necessary data
                name = compound_data.get('name', '')
                synonyms = compound_data.get('synonyms', [])
                inchi = compound_data.get('inchi', '')
                inchi_key = compound_data.get('inchi_key', '')
                smiles = compound_data.get('smiles', '')
                molecular_formula = compound_data.get('molecular_formula', '')
                molecular_weight = compound_data.get('molecular_weight', 0.0)
                
                # Convert properties to JSON
                properties = {
                    'pubchem': {
                        'logs': {
                            'imported_at': datetime.now().isoformat(),
                            'source': 'PubChem',
                            'version': '2.0'
                        },
                        'basic': {
                            'name': name,
                            'synonyms': synonyms[:10],  # Limit to 10 synonyms
                            'formula': molecular_formula,
                            'molecular_weight': molecular_weight
                        },
                        'identifiers': {
                            'cid': cid,
                            'inchi': inchi,
                            'inchi_key': inchi_key,
                            'smiles': smiles
                        },
                        'properties': self._extract_properties(compound_data)
                    }
                }
                
                with ConnectionManager() as conn:
                    with conn.cursor() as cursor:
                        if is_new:
                            # Insert new molecule
                            insert_query = """
                            INSERT INTO molecules 
                            (id, name, properties, created_at, updated_at) 
                            VALUES (%s, %s, %s, NOW(), NOW())
                            """
                            cursor.execute(insert_query, 
                                          (internal_id, name, json.dumps(properties)))
                            
                            logger.debug("Inserted new molecule %s (CID %s)", internal_id, cid)
                        else:
                            # Update existing molecule
                            update_query = """
                            UPDATE molecules 
                            SET name = %s, properties = properties || %s, updated_at = NOW() 
                            WHERE id = %s
                            """
                            cursor.execute(update_query, 
                                          (name, json.dumps(properties), internal_id))
                            
                            logger.debug("Updated molecule %s (CID %s)", internal_id, cid)
                        
                        # Commit transaction
                        conn.commit()
                
                # Successfully processed
                return
                
            except Exception as e:
                retry_count += 1
                conn.rollback()
                logger.warning("Database operation failed for CID %s (attempt %d/%d): %s", 
                              cid, retry_count, self.max_retries, str(e))
                time.sleep(1)
        
        # Failed after max retries
        logger.error("Failed to insert/update CID %s after %d attempts", cid, self.max_retries)
        raise Exception(f"Database operation failed for CID {cid}")
    
    def _extract_properties(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract useful properties from compound data.
        
        Args:
            compound_data: PubChem compound data
            
        Returns:
            Dictionary of extracted properties
        """
        properties = {}
        
        # Extract physical properties
        if 'props' in compound_data:
            for prop in compound_data.get('props', []):
                if prop.get('urn', {}).get('label') == 'LogP':
                    properties['logP'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Water Solubility':
                    properties['water_solubility'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Melting Point':
                    properties['melting_point'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Boiling Point':
                    properties['boiling_point'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Complexity':
                    properties['complexity'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'H-Bond Donor':
                    properties['h_bond_donor_count'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'H-Bond Acceptor':
                    properties['h_bond_acceptor_count'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'Rotatable Bond':
                    properties['rotatable_bond_count'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'Heavy Atom':
                    properties['heavy_atom_count'] = prop.get('value', {}).get('ival')
        
        return properties
    
    def _process_batch(self, cids: List[str]) -> None:
        """
        Process a batch of compounds in parallel.
        
        Args:
            cids: List of PubChem CIDs to process
        """
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.batch_size) as executor:
            # Submit all tasks
            future_to_cid = {executor.submit(self._process_compound, cid): cid for cid in cids}
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(future_to_cid):
                cid = future_to_cid[future]
                try:
                    success = future.result()
                    if success:
                        logger.debug("Successfully processed CID %s", cid)
                    else:
                        logger.warning("Failed to process CID %s", cid)
                except Exception as e:
                    logger.error("Exception processing CID %s: %s", cid, str(e))
                    self.checkpoint['failed_cids'][cid] = str(e)
                    self.stats['compounds_failed'] += 1
    
    def import_compounds(self, cid_list: List[str]) -> None:
        """
        Import compounds from a list of CIDs.
        
        Args:
            cid_list: List of PubChem CIDs to import
        """
        try:
            if not cid_list:
                logger.warning("Empty CID list provided, nothing to import")
                return
            
            # Update checkpoint with total compounds
            self.checkpoint['total_compounds'] = len(cid_list)
            
            # Start timing
            self.stats['start_time'] = time.time()
            
            # Process in batches
            total_batches = (len(cid_list) + self.batch_size - 1) // self.batch_size
            
            for batch_idx in range(total_batches):
                if batch_idx <= self.checkpoint.get('last_completed_batch', -1):
                    logger.info("Skipping batch %d/%d (already completed)", 
                                batch_idx + 1, total_batches)
                    continue
                
                # Get batch of CIDs
                start_idx = batch_idx * self.batch_size
                end_idx = min((batch_idx + 1) * self.batch_size, len(cid_list))
                batch_cids = cid_list[start_idx:end_idx]
                
                logger.info("Processing batch %d/%d (%d compounds)", 
                           batch_idx + 1, total_batches, len(batch_cids))
                
                # Process batch
                self._process_batch(batch_cids)
                
                # Update checkpoint
                self.checkpoint['last_completed_batch'] = batch_idx
                self._save_checkpoint()
                
                # Log progress
                progress = min(100.0, 100.0 * (batch_idx + 1) / total_batches)
                logger.info("Progress: %.1f%% (%d/%d batches, %d/%d compounds)", 
                           progress, batch_idx + 1, total_batches, 
                           len(self.checkpoint['processed_cids']), 
                           self.checkpoint['total_compounds'])
                
                # Sleep to avoid overtaxing PubChem API
                time.sleep(1)
            
            # Finish timing
            self.stats['end_time'] = time.time()
            self.stats['total_elapsed_time'] = self.stats['end_time'] - self.stats['start_time']
            
            # Final report
            logger.info("Import completed!")
            logger.info("Total compounds processed: %d", self.stats['compounds_processed'])
            logger.info("  - Added: %d", self.stats['compounds_added'])
            logger.info("  - Updated: %d", self.stats['compounds_updated'])
            logger.info("  - Skipped: %d", self.stats['compounds_skipped'])
            logger.info("  - Failed: %d", self.stats['compounds_failed'])
            logger.info("Total elapsed time: %.2f seconds", self.stats['total_elapsed_time'])
            
            # Save final checkpoint
            self._save_checkpoint()
            
        except Exception as e:
            logger.error("Import failed: %s", str(e))
            # Save checkpoint
            self._save_checkpoint()
    
    def generate_report(self, output_file: str) -> None:
        """
        Generate import report.
        
        Args:
            output_file: Path to output report file
        """
        try:
            # Create report data
            report = {
                'timestamp': datetime.now().isoformat(),
                'statistics': self.stats,
                'processed_compounds': len(self.checkpoint['processed_cids']),
                'failed_compounds': len(self.checkpoint['failed_cids']),
                'total_compounds': self.checkpoint['total_compounds'],
                'completion_percentage': (len(self.checkpoint['processed_cids']) / 
                                         max(1, self.checkpoint['total_compounds']) * 100),
                'failed_compound_details': self.checkpoint['failed_cids']
            }
            
            # Calculate rates
            if self.stats['total_elapsed_time'] > 0:
                report['processing_rate'] = (self.stats['compounds_processed'] / 
                                          self.stats['total_elapsed_time'])
            else:
                report['processing_rate'] = 0
                
            # Save report
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)
            
            logger.info("Report saved to %s", output_file)
            
        except Exception as e:
            logger.error("Failed to generate report: %s", str(e))

def main():
    """Main function for PubChem import."""
    parser = argparse.ArgumentParser(description='Import PubChem compounds to Supabase')
    parser.add_argument('--cid-file', type=str, help='File containing PubChem CIDs (one per line)')
    parser.add_argument('--cid-list', type=str, help='Comma-separated list of PubChem CIDs')
    parser.add_argument('--batch-size', type=int, default=BATCH_SIZE, 
                       help='Number of compounds to process in parallel')
    parser.add_argument('--checkpoint-file', type=str, default=CHECKPOINT_FILE,
                       help='Path to checkpoint file')
    parser.add_argument('--report-file', type=str, 
                       default=f"reports/pubchem_import_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                       help='Path to output report file')
    args = parser.parse_args()
    
    # Get CID list
    cid_list = []
    if args.cid_file:
        try:
            with open(args.cid_file, 'r') as f:
                cid_list = [line.strip() for line in f if line.strip()]
            logger.info("Loaded %d CIDs from %s", len(cid_list), args.cid_file)
        except Exception as e:
            logger.error("Failed to load CID file: %s", str(e))
            return 1
    elif args.cid_list:
        cid_list = [cid.strip() for cid in args.cid_list.split(',') if cid.strip()]
        logger.info("Using %d CIDs from command line", len(cid_list))
    else:
        # Use CIDs from identifier manager
        id_manager = CryoprotectantIdentifierManager.get_instance()
        cid_list = id_manager.get_all_pubchem_cids()
        logger.info("Using %d CIDs from identifier manager", len(cid_list))
    
    if not cid_list:
        logger.error("No CIDs provided. Use --cid-file or --cid-list option.")
        return 1
    
    # Create importer
    importer = PubChemImporter(
        checkpoint_file=args.checkpoint_file,
        batch_size=args.batch_size
    )
    
    # Import compounds
    importer.import_compounds(cid_list)
    
    # Generate report
    importer.generate_report(args.report_file)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

#### Step 3.2: Create Test Script for Enhanced PubChem Import

Create a simple test script to verify the enhanced importer works correctly:

```python
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
```

#### Step 3.3: Run Test Import and Verify Results

```bash
# Run test import
python test_pubchem_enhanced_import.py

# Verify database entries
# Check that molecules were inserted correctly
```

### Success Criteria
- √ Parallel batch processing of PubChem compounds
- √ Automatic resumption from checkpoint if interrupted
- √ Integration with connection pool for efficient DB connections
- √ Integration with identifier manager for consistent molecule IDs
- √ Comprehensive error handling with retry logic
- √ Detailed statistics and reporting
- √ Performance improvement of at least 40% compared to original script

## 4. ChEMBL Import Enhancement (Task 2.7)

This follows a similar pattern to the PubChem Import Optimization, with specific adaptations for the ChEMBL API and data model. The implementation details would be documented here.

## 5. Molecule Viewer Optimization (Task 3.1)

This involves enhancing the molecule viewer component in the web interface to efficiently display the imported scientific data, with optimizations for caching and rendering.

## 6. Data Table Enhancement (Task 3.3)

This involves improving the data tables in the web interface to display comprehensive molecule properties with sorting, filtering, and pagination.

## Implementation Sequence

### Week 1
1. Implement Connection Pool Enhancement (Task 1.1)
2. Create Cryoprotectant Identifier List (Task 2.2)
3. Implement PubChem Import Optimization (Task 2.3)

### Week 2
4. Implement ChEMBL Import Enhancement (Task 2.7)
5. Optimize Molecule Viewer (Task 3.1)
6. Enhance Data Tables (Task 3.3)

## Progress Tracking

Track implementation progress using the following format:

| Task ID | Description | Status | Completion Date | Notes |
|---------|-------------|--------|----------------|-------|
| 1.1 | Connection Pool Enhancement | Not Started | | |
| 2.2 | Cryoprotectant Identifier List | Not Started | | |
| 2.3 | PubChem Import Optimization | Not Started | | |
| 2.7 | ChEMBL Import Enhancement | Not Started | | |
| 3.1 | Molecule Viewer Optimization | Not Started | | |
| 3.3 | Data Table Enhancement | Not Started | | |