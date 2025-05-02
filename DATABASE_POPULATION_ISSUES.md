# Database Population Issues & Optimization Guide

## Executive Summary

Our database population task has encountered several persistent issues that need to be addressed before proceeding with large-scale data import. This document identifies key problems, root causes, and provides concrete solutions for Roo Code agents to implement.

Based on our testing, the Supabase session pooler connection (port 6543) has proven to be the most reliable connection method and should be prioritized in our implementation strategy. The document has been updated to reflect this finding while maintaining fallback mechanisms for maximum resilience.

## 1. Critical Connection Issues

### 1.1 Problem: Connection Timeouts and Instability
- Connection attempts to Supabase frequently timeout
- DNS resolution issues between hostname and IP address
- Transaction management issues leading to "current transaction is aborted" errors

### 1.2 Root Causes
- Supabase DNS resolution is unreliable in certain network environments
- Connection pooling not properly handling dropped connections
- Missing transaction rollback in error cases causing cascading failures

### 1.3 Solutions for Implementation
```python
# 1. Use session pooler connection as primary method (most reliable)
# with explicit IP-based connection as fallback
def get_db_connection():
    # First try session pooler connection (most reliable)
    try:
        db_host = os.getenv('SUPABASE_DB_HOST')
        db_port = os.getenv('SUPABASE_DB_PORT', '6543')  # Session pooler port
        db_name = os.getenv('SUPABASE_DB_NAME')
        db_user = os.getenv('SUPABASE_DB_USER')
        db_password = os.getenv('SUPABASE_DB_PASSWORD')
        
        conn = psycopg2.connect(
            host=db_host,
            port=db_port,
            dbname=db_name,
            user=db_user,
            password=db_password
        )
        logger.info(f"Connected to session pooler on {db_host}:{db_port}")
        return conn
    except Exception as e:
        logger.warning(f"Session pooler connection failed: {e}")
    
    # Next try IP-based direct connection
    db_ip = os.getenv('SUPABASE_DB_IP_ADDRESS')
    if db_ip:
        try:
            return connect_with_ip(db_ip)
        except Exception as e:
            logger.warning(f"IP-based connection failed: {e}")
    
    # Fall back to hostname with resolution
    try:
        db_host = os.getenv('SUPABASE_DB_HOST')
        ip = resolve_hostname(db_host)
        return connect_with_ip(ip)
    except Exception:
        # Last resort: try MCP-based connection
        return get_mcp_connection()
        
# 2. Implement circuit breaker pattern for connection stability
class CircuitBreaker:
    def __init__(self, failure_threshold=3, recovery_timeout=60):
        self.failure_count = 0
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.last_failure_time = 0
        self.state = "CLOSED"  # CLOSED, OPEN, HALF-OPEN
        
    def execute(self, func, *args, **kwargs):
        if self.state == "OPEN":
            if time.time() - self.last_failure_time > self.recovery_timeout:
                self.state = "HALF-OPEN"
            else:
                raise ConnectionError("Circuit breaker open")
                
        try:
            result = func(*args, **kwargs)
            if self.state == "HALF-OPEN":
                self.state = "CLOSED"
                self.failure_count = 0
            return result
        except Exception as e:
            self.failure_count += 1
            self.last_failure_time = time.time()
            if self.failure_count >= self.failure_threshold:
                self.state = "OPEN"
            raise
```

## 2. Transaction Management Issues

### 2.1 Problem: Transaction Leaks and Deadlocks
- Incomplete transaction handling leading to connection pool exhaustion
- Nested transaction attempts causing errors
- Deadlocks during concurrent property operations

### 2.2 Root Causes
- Missing transaction cleanup in error cases
- No detection of existing transactions before starting new ones
- Concurrent operations on the same property types

### 2.3 Solutions for Implementation
```python
# 1. Safe transaction context manager with cleanup
@contextmanager
def safe_transaction():
    connection = get_db_connection()
    transaction = None
    try:
        # Check if already in transaction
        status_query = "SELECT current_setting('transaction_isolation')"
        try:
            connection.execute_query(status_query)
            # If we get here, we're not in an aborted transaction
        except Exception as e:
            if "current transaction is aborted" in str(e):
                logger.warning("Detected aborted transaction, rolling back")
                connection.rollback()
        
        # Start new transaction
        transaction = connection.begin_transaction()
        yield transaction
        connection.commit_transaction(transaction)
    except Exception as e:
        if transaction:
            try:
                connection.rollback_transaction(transaction)
            except Exception as rollback_error:
                logger.error(f"Rollback failed: {rollback_error}")
        raise
    finally:
        # Verify transaction cleanup
        try:
            connection.execute_query("SELECT 1")
        except Exception:
            # If query fails, attempt emergency rollback
            try:
                connection.rollback()
            except:
                pass

# 2. Lock-based property type creation to prevent race conditions
property_type_locks = {}
def get_or_create_property_type(name, data_type='numeric'):
    # Use a lock per property name to prevent race conditions
    if name not in property_type_locks:
        property_type_locks[name] = threading.Lock()
    
    with property_type_locks[name]:
        # Rest of implementation...
```

## 3. Batch Processing Inefficiencies

### 3.1 Problem: Inefficient Batch Processing
- Individual property inserts instead of bulk operations
- Batch sizes too large causing timeouts or memory issues
- No checkpointing for resumable operations

### 3.2 Root Causes
- Transaction overhead for each property
- Missing bulk insert optimizations
- No state persistence between runs

### 3.3 Solutions for Implementation
```python
# 1. Optimized property batch insertion
def bulk_insert_properties(property_records, batch_size=1000):
    """Insert multiple property records in optimized batches"""
    if not property_records:
        return 0
        
    # Group by property type for more efficient inserts
    grouped_properties = {}
    for record in property_records:
        prop_type = record['property_type_id']
        if prop_type not in grouped_properties:
            grouped_properties[prop_type] = []
        grouped_properties[prop_type].append(record)
    
    # Insert each property type group in a single transaction
    total_inserted = 0
    for prop_type, records in grouped_properties.items():
        # Process in reasonable batch sizes
        for i in range(0, len(records), batch_size):
            batch = records[i:i+batch_size]
            with safe_transaction():
                # Build multi-value insert for this batch
                columns = ['molecule_id', 'property_type_id', 'numeric_value', 
                           'text_value', 'boolean_value', 'created_by']
                values_list = []
                params = []
                
                for record in batch:
                    values_list.append(f"(%s, %s, %s, %s, %s, %s)")
                    params.extend([
                        record['molecule_id'],
                        record['property_type_id'],
                        record.get('numeric_value'),
                        record.get('text_value'),
                        record.get('boolean_value'),
                        record.get('created_by')
                    ])
                
                # Execute the multi-value insert
                query = f"""
                INSERT INTO molecular_properties 
                ({', '.join(columns)})
                VALUES {', '.join(values_list)}
                """
                execute_query(query, params)
                total_inserted += len(batch)
    
    return total_inserted

# 2. Resumable import with checkpointing
def resumable_batch_import(items, process_func, checkpoint_file, batch_size=100):
    """Process items in batches with checkpoint-based resumability"""
    # Load checkpoint if exists
    checkpoint = {}
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
        except:
            logger.warning(f"Could not load checkpoint from {checkpoint_file}")
    
    # Get starting position
    position = checkpoint.get('position', 0)
    processed_count = checkpoint.get('processed', 0)
    
    # Process remaining items
    total_batches = math.ceil((len(items) - position) / batch_size)
    start_time = time.time()
    
    for batch_num in range(total_batches):
        batch_start = position + (batch_num * batch_size)
        batch_end = min(batch_start + batch_size, len(items))
        batch = items[batch_start:batch_end]
        
        # Process this batch
        logger.info(f"Processing batch {batch_num+1}/{total_batches}, items {batch_start}-{batch_end}")
        batch_result = process_func(batch)
        
        # Update checkpoint
        processed_count += len(batch)
        checkpoint = {
            'position': batch_end,
            'processed': processed_count,
            'last_updated': datetime.now().isoformat(),
            'last_batch': batch_num
        }
        
        # Save checkpoint
        with open(checkpoint_file, 'w') as f:
            json.dump(checkpoint, f)
            
        # Report progress
        elapsed = time.time() - start_time
        items_per_sec = processed_count / elapsed if elapsed > 0 else 0
        logger.info(f"Progress: {processed_count}/{len(items)} items "
                   f"({processed_count/len(items)*100:.1f}%) "
                   f"at {items_per_sec:.1f} items/sec")
    
    return processed_count
```

## 4. ChEMBL Data Import Strategy Issues

### 4.1 Problem: Search Term Approach Not Finding Compounds
- Search terms like "cryoprotectant" returning 0 compounds
- Reference compounds list not properly utilized
- Missing property-based filtering

### 4.2 Root Causes
- ChEMBL doesn't categorize compounds by application
- Search focused on names/descriptions rather than chemical properties
- Missing similarity search capabilities

### 4.3 Solutions for Implementation
```python
# 1. Property-based cryoprotectant identification
def find_potential_cryoprotectants():
    """Find compounds with properties matching typical cryoprotectants"""
    # Cryoprotectants typically have:
    # - Multiple hydrogen bond donors/acceptors
    # - Reasonable water solubility (LogP)
    # - Molecular weight in specific range
    query = """
    SELECT molecule_chembl_id, canonical_smiles, pref_name, molecule_properties
    FROM molecule_dictionary md
    JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
    JOIN compound_properties cp ON md.molregno = cp.molregno
    WHERE 
        cp.mw_freebase BETWEEN 30 AND 500
        AND cp.alogp BETWEEN -3 AND 3
        AND cp.hba BETWEEN 2 AND 20
        AND cp.hbd BETWEEN 1 AND 10
        AND md.max_phase > 0
    ORDER BY md.max_phase DESC
    LIMIT 5000
    """
    return chembl_client.execute_query(query)

# 2. Structural similarity search from known cryoprotectants
def find_similar_compounds(reference_smiles, similarity_threshold=70):
    """Find compounds similar to known cryoprotectants"""
    similar_molecules = []
    for smiles in reference_smiles:
        results = chembl_client.molecule.filter(
            similarity=smiles,
            similarity_threshold=similarity_threshold
        ).only(['molecule_chembl_id', 'molecule_structures'])
        similar_molecules.extend(results)
    return similar_molecules

# 3. Chemical class filtering using SMARTS patterns
def identify_compounds_by_chemical_class():
    """Find compounds belonging to chemical classes common for cryoprotectants"""
    # SMARTS patterns for different cryoprotectant classes
    chemical_classes = {
        'polyols': '[OX2H][CX4][CX4][OX2H]', # Pattern for compounds with adjacent hydroxyl groups
        'amides': '[NX3][CX3]=[OX1]',        # Pattern for amide functional group
        'sulfoxides': '[#16X3]=[OX1]',       # Pattern for sulfoxide group (like DMSO)
        'sugars': '[OX2H][CX4][CX4][CX4][OX2H]' # Simplified sugar pattern
    }
    
    results = {}
    for class_name, smarts in chemical_classes.items():
        query = f"""
        SELECT molecule_chembl_id, canonical_smiles, pref_name
        FROM molecule_dictionary md
        JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
        JOIN compound_structures cs ON mh.parent_molregno = cs.molregno
        WHERE mol_to_smarts(cs.molfile::mol) @> '{smarts}'::smarts
        LIMIT 1000
        """
        results[class_name] = chembl_client.execute_query(query)
    
    return results
```

## 5. Integration & Strategic Recommendations

### 5.1 Recommended Overall Approach
1. **Implement Connection Resilience Layer**
   - Prioritize session pooler connection (most reliable approach per testing)
   - Fall back to direct IP connections and MCP methods when needed
   - Add circuit breaker pattern to prevent cascading failures
   - Implement robust transaction handling with cleanup guarantees

2. **Optimize Batch Processing**
   - Add checkpointing for resumability
   - Implement grouped property inserts
   - Use dynamically adjusted batch sizes based on performance

3. **Enhance ChEMBL Data Import**
   - Switch to property-based filtering approach
   - Implement similarity search from known reference compounds
   - Add chemical class-based filtering with SMARTS patterns

### 5.2 Implementation Strategy
Use a phased implementation approach:
1. **Phase 1**: Fix connection and transaction handling issues
2. **Phase 2**: Implement optimized batch processing
3. **Phase 3**: Enhance ChEMBL data discovery methods
4. **Phase 4**: Integrate all components into a resilient population system

### 5.3 Monitoring & Verification
After implementation, establish:
1. Connection health monitoring with automatic fallback
2. Progress tracking with time estimation
3. Data quality verification at checkpoints
4. Performance metrics collection

## Conclusion

The current database population issues stem primarily from connection instability, transaction management problems, and inefficient batch processing. By implementing the solutions outlined in this document, we can create a resilient, efficient, and reliable data import system that overcomes these challenges.

Our testing confirms that the Supabase session pooler connection (port 6543) offers the best stability and performance for database operations. The implementation strategy prioritizes this connection method while maintaining a robust fallback system for maximum resilience.

For the Roo Code agents, this document provides the specific implementation details needed to fix these issues and move forward with successfully populating the database with both PubChem and ChEMBL data.