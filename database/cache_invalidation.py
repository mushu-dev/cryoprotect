#!/usr/bin/env python3
"""
Cache invalidation strategies for database updates.

This module provides functions and hooks to invalidate the cache when data
is modified in the database. It uses a combination of triggers, database
event listeners, and explicit invalidation calls to keep the cache synchronized
with the database.
"""

import logging
from typing import Dict, List, Any, Optional, Set

# Import caching module
from database.cache import (
    invalidate_molecule,
    invalidate_property,
    cache_invalidate_pattern,
    clear_cache
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Tables that should trigger cache invalidation
CACHE_SENSITIVE_TABLES = {
    'molecules': {
        'pattern': 'molecule:*',
        'related_patterns': ['query:*']
    },
    'molecular_properties': {
        'pattern': 'property:*',
        'related_patterns': ['molecule:*', 'query:*']
    },
    'property_types': {
        'pattern': 'property:*',
        'related_patterns': ['query:*']
    },
    'property_calculation_queue': {
        'pattern': 'calculation:*',
        'related_patterns': []
    },
    'cryoprotection_scores': {
        'pattern': 'score:*',
        'related_patterns': ['molecule:*', 'query:*']
    },
    'mixtures': {
        'pattern': 'mixture:*',
        'related_patterns': ['query:*']
    }
}

# Map of dependency relationships between tables
TABLE_DEPENDENCIES = {
    'molecules': ['molecular_properties', 'cryoprotection_scores', 'mixtures'],
    'molecular_properties': ['cryoprotection_scores'],
    'property_types': ['molecular_properties'],
    'mixtures': []
}

def invalidate_table_cache(table_name: str) -> int:
    """
    Invalidate cache for a specific table.
    
    Args:
        table_name: The table name
        
    Returns:
        Number of keys invalidated
    """
    if table_name not in CACHE_SENSITIVE_TABLES:
        logger.debug(f"Table {table_name} is not cache sensitive, skipping invalidation")
        return 0
    
    table_info = CACHE_SENSITIVE_TABLES[table_name]
    pattern = table_info['pattern']
    related_patterns = table_info['related_patterns']
    
    count = cache_invalidate_pattern(pattern)
    
    # Also invalidate related patterns
    for related_pattern in related_patterns:
        count += cache_invalidate_pattern(related_pattern)
    
    logger.info(f"Invalidated {count} cache keys for table {table_name}")
    return count

def invalidate_record_cache(table_name: str, record_id: str) -> int:
    """
    Invalidate cache for a specific record.
    
    Args:
        table_name: The table name
        record_id: The record ID
        
    Returns:
        Number of keys invalidated
    """
    if table_name not in CACHE_SENSITIVE_TABLES:
        logger.debug(f"Table {table_name} is not cache sensitive, skipping invalidation")
        return 0
    
    count = 0
    
    # Special handling for specific tables
    if table_name == 'molecules':
        count = invalidate_molecule(record_id)
    elif table_name == 'molecular_properties':
        # For properties, we need the molecule_id as well
        try:
            from database.db import execute_query
            result = execute_query(
                "SELECT molecule_id FROM molecular_properties WHERE id = %s",
                (record_id,)
            )
            if result and 'molecule_id' in result[0]:
                molecule_id = result[0]['molecule_id']
                count = invalidate_molecule(molecule_id)
                count += invalidate_property()
        except Exception as e:
            logger.warning(f"Error getting molecule_id for property {record_id}: {e}")
            # Fallback to invalidating all property cache
            count = invalidate_property()
    else:
        # Generic approach for other tables
        table_info = CACHE_SENSITIVE_TABLES[table_name]
        pattern = table_info['pattern'].replace('*', f'*{record_id}*')
        count = cache_invalidate_pattern(pattern)
    
    logger.info(f"Invalidated {count} cache keys for {table_name} record {record_id}")
    return count

def create_invalidation_sql() -> str:
    """
    Create SQL statements for database triggers that invalidate cache.
    
    Returns:
        SQL string with trigger definitions
    """
    # We can't directly invalidate the cache from a PostgreSQL trigger,
    # but we can insert a record into a special table that our application
    # will monitor for invalidation events.
    
    sql = """
    -- Create cache invalidation table
    CREATE TABLE IF NOT EXISTS cache_invalidation_events (
        id SERIAL PRIMARY KEY,
        table_name VARCHAR(255) NOT NULL,
        record_id VARCHAR(255),
        operation VARCHAR(50) NOT NULL,
        created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
        processed_at TIMESTAMP WITH TIME ZONE DEFAULT NULL
    );
    
    -- Create function to record invalidation events
    CREATE OR REPLACE FUNCTION record_cache_invalidation()
    RETURNS TRIGGER AS $$
    BEGIN
        IF (TG_OP = 'INSERT') THEN
            INSERT INTO cache_invalidation_events 
                (table_name, record_id, operation) 
            VALUES 
                (TG_TABLE_NAME, NEW.id, 'INSERT');
            RETURN NEW;
        ELSIF (TG_OP = 'UPDATE') THEN
            INSERT INTO cache_invalidation_events 
                (table_name, record_id, operation) 
            VALUES 
                (TG_TABLE_NAME, NEW.id, 'UPDATE');
            RETURN NEW;
        ELSIF (TG_OP = 'DELETE') THEN
            INSERT INTO cache_invalidation_events 
                (table_name, record_id, operation) 
            VALUES 
                (TG_TABLE_NAME, OLD.id, 'DELETE');
            RETURN OLD;
        END IF;
        RETURN NULL;
    END;
    $$ LANGUAGE plpgsql;
    """
    
    # Create triggers for each cache-sensitive table
    for table_name in CACHE_SENSITIVE_TABLES:
        sql += f"""
        -- Trigger for {table_name}
        DROP TRIGGER IF EXISTS {table_name}_cache_invalidation ON {table_name};
        CREATE TRIGGER {table_name}_cache_invalidation
        AFTER INSERT OR UPDATE OR DELETE ON {table_name}
        FOR EACH ROW
        EXECUTE FUNCTION record_cache_invalidation();
        """
    
    return sql

class CacheInvalidationProcessor:
    """
    Class for processing cache invalidation events from the database.
    """
    def __init__(self):
        """Initialize the processor."""
        self.last_processed_id = 0
        self.invalidated_tables = set()
    
    def process_events(self, batch_size: int = 100) -> int:
        """
        Process cache invalidation events from the database.
        
        Args:
            batch_size: Maximum number of events to process in a batch
            
        Returns:
            Number of events processed
        """
        try:
            from database.db import execute_query, transaction
            
            # Get pending events
            query = """
                SELECT id, table_name, record_id, operation 
                FROM cache_invalidation_events 
                WHERE processed_at IS NULL AND id > %s
                ORDER BY id 
                LIMIT %s
            """
            events = execute_query(query, (self.last_processed_id, batch_size))
            
            if not events:
                return 0
            
            # Process each event
            processed_count = 0
            processed_ids = []
            
            for event in events:
                event_id = event['id']
                table_name = event['table_name']
                record_id = event['record_id']
                operation = event['operation']
                
                try:
                    # Invalidate cache for this event
                    if record_id:
                        invalidate_record_cache(table_name, record_id)
                    else:
                        invalidate_table_cache(table_name)
                    
                    # Keep track of invalidated tables
                    self.invalidated_tables.add(table_name)
                    
                    # Mark event as processed
                    processed_ids.append(event_id)
                    processed_count += 1
                    
                    # Update last processed ID
                    if event_id > self.last_processed_id:
                        self.last_processed_id = event_id
                        
                except Exception as e:
                    logger.error(f"Error processing invalidation event {event_id}: {e}")
            
            # Mark events as processed in the database
            if processed_ids:
                with transaction() as cursor:
                    cursor.execute(
                        "UPDATE cache_invalidation_events SET processed_at = NOW() WHERE id = ANY(%s)",
                        (processed_ids,)
                    )
            
            logger.info(f"Processed {processed_count} cache invalidation events")
            return processed_count
            
        except Exception as e:
            logger.error(f"Error processing cache invalidation events: {e}")
            return 0
    
    def invalidate_dependent_tables(self) -> int:
        """
        Invalidate cache for tables dependent on modified tables.
        
        Returns:
            Number of tables invalidated
        """
        if not self.invalidated_tables:
            return 0
        
        dependent_tables = set()
        for table in self.invalidated_tables:
            if table in TABLE_DEPENDENCIES:
                dependent_tables.update(TABLE_DEPENDENCIES[table])
        
        # Remove tables already invalidated
        dependent_tables -= self.invalidated_tables
        
        # Invalidate dependent tables
        count = 0
        for table in dependent_tables:
            count += invalidate_table_cache(table)
        
        # Reset invalidated tables
        self.invalidated_tables = set()
        
        return len(dependent_tables)
    
    def cleanup_old_events(self, days: int = 7) -> int:
        """
        Clean up old invalidation events.
        
        Args:
            days: Delete events older than this many days
            
        Returns:
            Number of events deleted
        """
        try:
            from database.db import execute_query
            
            query = """
                DELETE FROM cache_invalidation_events 
                WHERE processed_at IS NOT NULL 
                  AND created_at < NOW() - INTERVAL '%s days'
                RETURNING COUNT(*) as count
            """
            result = execute_query(query, (days,))
            
            if result and 'count' in result[0]:
                count = result[0]['count']
                logger.info(f"Deleted {count} old cache invalidation events")
                return count
            
            return 0
            
        except Exception as e:
            logger.error(f"Error cleaning up old invalidation events: {e}")
            return 0

# Create a singleton instance of the processor
invalidation_processor = CacheInvalidationProcessor()

def setup_cache_invalidation() -> bool:
    """
    Set up cache invalidation infrastructure.
    
    Returns:
        True if successful, False otherwise
    """
    try:
        from database.db import execute_query
        
        # Create the invalidation table and triggers
        sql = create_invalidation_sql()
        execute_query(sql)
        
        logger.info("Cache invalidation infrastructure created successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error setting up cache invalidation: {e}")
        return False

def process_invalidation_events() -> int:
    """
    Process pending cache invalidation events.
    
    Returns:
        Number of events processed
    """
    count = invalidation_processor.process_events()
    if count > 0:
        invalidation_processor.invalidate_dependent_tables()
    return count