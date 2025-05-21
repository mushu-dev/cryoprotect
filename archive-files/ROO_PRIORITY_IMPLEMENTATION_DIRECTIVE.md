# ROO PRIORITY IMPLEMENTATION DIRECTIVE

## TASK OVERVIEW
Implement the first phase of data flow optimization and scientific data integration for CryoProtect v2. This phase focuses on three critical components:

1. Connection Pool Enhancement
2. Cryoprotectant Identifier List
3. PubChem Import Optimization

## SUCCESS CRITERIA
- Connection pooling demonstrates 30%+ performance improvement in data import
- Master identifier list contains at least 50 common cryoprotectants with complete identifiers
- PubChem import process shows 40%+ performance improvement with parallel processing
- All implementations have comprehensive error handling and logging
- Data imports are resumable from checkpoints if interrupted

## FILE REFERENCES

### Key Files to Create:
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/connection_pool_wrapper.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/cryoprotectant_identifiers.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/data/cryoprotectant_master_list.json`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/PubChem_CryoProtectants_Supabase_Enhanced.py`

### Key Files to Modify:
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/config.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/pubchem/client.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/pubchem/cache.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/pubchem/rate_limiter.py`

## IMPLEMENTATION TASKS

### Task 1: Connection Pool Enhancement
Create a robust connection pooling system to reduce database connection overhead and improve data import performance.

**Specific Implementation Instructions:**

1. Create `connection_pool_wrapper.py` with the following components:
   - ConnectionPoolWrapper class with singleton pattern
   - Automatic health checking and reconnection
   - Connection lifecycle management
   - Context manager for safe connection handling

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

2. Update `config.py` to support connection pooling:

```python
# Add to config.py
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

3. Create a test script to verify connection pooling:

```python
# test_connection_pool.py
import logging
import sys
import threading
import time
from connection_pool_wrapper import ConnectionPoolWrapper, ConnectionManager
from config import get_db_config

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run_test_query(thread_id):
    """Run a test query in a thread."""
    try:
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                cursor.execute("SELECT pg_sleep(0.5)")
                logger.info(f"Thread {thread_id} completed query")
    except Exception as e:
        logger.error(f"Thread {thread_id} query failed: {str(e)}")

def test_connection_pool():
    """Test the connection pool with concurrent connections."""
    # Create connection pool
    config = get_db_config()
    pool = ConnectionPoolWrapper.get_instance(config)
    
    # Test basic query
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            cursor.execute("SELECT 1 as test")
            result = cursor.fetchone()
            logger.info(f"Basic query test: {result}")
    
    # Test multiple connections
    threads = []
    for i in range(5):
        thread = threading.Thread(target=run_test_query, args=(i,))
        threads.append(thread)
        thread.start()
    
    # Wait for threads to complete
    for thread in threads:
        thread.join()
    
    # Check active connections
    logger.info(f"Active connections: {pool.active_connections}")
    
    # Test health check
    pool._check_pool_health()
    
    return True

if __name__ == "__main__":
    logger.info("Testing connection pool...")
    test_connection_pool()
    logger.info("Connection pool test completed!")
```

### Task 2: Cryoprotectant Identifier List
Create a standardized list of cryoprotectant molecule identifiers to ensure consistent data import across multiple data sources.

**Specific Implementation Instructions:**

1. Create the `cryoprotectant_identifiers.py` module:

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

2. Create a test script for the identifier manager:

```python
# test_cryoprotectant_identifiers.py
import logging
from cryoprotectant_identifiers import CryoprotectantIdentifierManager, initialize_cryoprotectant_list

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_identifier_manager():
    """Test the cryoprotectant identifier manager."""
    # Initialize with common cryoprotectants
    initialize_cryoprotectant_list()
    
    # Get instance
    manager = CryoprotectantIdentifierManager.get_instance()
    
    # Test lookup by different identifiers
    glycerol_by_cid = manager.get_internal_id_by_pubchem_cid("962")
    glycerol_by_name = manager.get_internal_id_by_name("Glycerol")
    glycerol_by_inchi = manager.get_internal_id_by_inchi_key("PEDCQBHIVMGVHV-UHFFFAOYSA-N")
    
    logger.info(f"Glycerol by CID: {glycerol_by_cid}")
    logger.info(f"Glycerol by name: {glycerol_by_name}")
    logger.info(f"Glycerol by InChI key: {glycerol_by_inchi}")
    
    # Test resolve_identifier
    dmso_id, confidence = manager.resolve_identifier(pubchem_cid="5988")
    logger.info(f"DMSO resolved to {dmso_id} with confidence {confidence}")
    
    # Test adding a new molecule
    manager.add_molecule("CRYO006", {
        "internal_id": "CRYO006",
        "pubchem_cid": "702",
        "chembl_id": "CHEMBL388978",
        "cas_number": "50-70-4",
        "names": ["Sorbitol", "D-Glucitol", "D-Sorbitol"],
        "inchi_key": "FBPFZTCFMRRESA-JGWLITMVSA-N",
        "smiles": "C(C(C(C(C(CO)O)O)O)O)O",
        "formula": "C6H14O6",
        "molecular_weight": 182.17,
        "category": "polyol"
    })
    
    # Verify the new molecule
    sorbitol_id = manager.get_internal_id_by_name("Sorbitol")
    logger.info(f"Sorbitol ID: {sorbitol_id}")
    
    # Save the updated list
    manager.save_identifiers()
    
    return True

if __name__ == "__main__":
    logger.info("Testing cryoprotectant identifier manager...")
    test_identifier_manager()
    logger.info("Cryoprotectant identifier test completed!")
```

### Task 3: PubChem Import Optimization
Create an enhanced PubChem importer with parallel batch processing, robust error handling, and integration with the connection pool and identifier manager.

**Specific Implementation Instructions:**

1. Create `PubChem_CryoProtectants_Supabase_Enhanced.py`:

```python
# PubChem_CryoProtectants_Supabase_Enhanced.py
import os
import json
import logging
import time
import argparse
import sys
import concurrent.futures
import threading
from typing import Dict, List, Any, Optional, Set, Tuple
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

2. Create a test script for the enhanced PubChem importer:

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

## TESTING INSTRUCTIONS

### Testing Connection Pool

1. Run the connection pool test:
```bash
python test_connection_pool.py
```

2. Verify the following output:
   - Connection pool initialized successfully
   - Basic query test passes
   - Multiple concurrent connections work correctly
   - Health check passes

### Testing Cryoprotectant Identifier Manager

1. Run the identifier manager test:
```bash
python test_cryoprotectant_identifiers.py
```

2. Verify the following output:
   - Identifier manager initialized with common cryoprotectants
   - Lookups by different identifier types work correctly
   - New molecule can be added successfully

### Testing PubChem Enhanced Import

1. Run the PubChem import test with a small set of CIDs:
```bash
python test_pubchem_enhanced_import.py --cids 962,5988,6342
```

2. Verify the following output:
   - CIDs are processed successfully
   - Checkpointing works correctly
   - Database entries are created/updated
   - Performance report shows improvement over original implementation

## INTEGRATION SEQUENCE

1. First implement and test the connection pool implementation
2. Next implement and test the cryoprotectant identifier manager
3. Finally implement and test the enhanced PubChem importer

Each component builds on the previous one, so implement them in sequence.

## VERIFICATION STEPS

After implementation, verify the following:

1. Using the connection pool reduces connection overhead by at least 30%
2. The cryoprotectant identifier list contains all common molecules with complete identifiers
3. The enhanced PubChem importer shows at least 40% performance improvement
4. All implementations have comprehensive error handling and logging
5. Data imports are resumable from checkpoints if interrupted

Once verified, the implementation can proceed to the next phase.