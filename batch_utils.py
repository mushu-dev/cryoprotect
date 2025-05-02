#!/usr/bin/env python3
"""
Batch Processing Utilities for CryoProtect v2

This module provides optimized batch processing utilities for database operations:
1. Efficient property batch insertion with grouping by property type
2. Transaction-safe batch operations with proper error handling
3. Optimized multi-value inserts for improved performance

Based on specifications in DATABASE_POPULATION_ISSUES.md (Section 3.3.1)
"""

import logging
import math
import os
import json
import time
from typing import List, Dict, Any, Optional, Union, Tuple
from datetime import datetime
from contextlib import contextmanager
from uuid import UUID

# Import transaction_utils and sql_executor with error handling
try:
    from transaction_utils import safe_transaction
    from sql_executor import execute_query, process_in_batches
except ImportError:
    # Mock implementations for testing
    from contextlib import contextmanager
    
    @contextmanager
    def safe_transaction():
        """Mock implementation of safe_transaction for testing."""
        yield None
    
    def execute_query(query, params=None, fetch_one=False, dict_cursor=True):
        """Mock implementation of execute_query for testing."""
        return []
    
    def process_in_batches(items, batch_size=1000, process_func=None):
        """Mock implementation of process_in_batches for testing."""
        if process_func and items:
            for i in range(0, len(items), batch_size):
                process_func(items[i:i+batch_size])
        return []

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('batch_utils')


def bulk_insert_properties(property_records: List[Dict[str, Any]], batch_size: int = 1000) -> int:
    """
    Insert multiple property records in optimized batches.
    
    This function implements the optimized property batch insertion algorithm
    specified in DATABASE_POPULATION_ISSUES.md (Section 3.3.1). It groups
    property records by property type for more efficient inserts and processes
    them in reasonable batch sizes to avoid timeouts or memory issues.
    
    Args:
        property_records: List of property records to insert. Each record should be a dictionary
                         with the following keys:
                         - molecule_id: UUID of the molecule
                         - property_type_id: UUID of the property type
                         - numeric_value: Optional numeric value
                         - text_value: Optional text value
                         - boolean_value: Optional boolean value
                         - created_by: Optional UUID of the user creating the property
        batch_size: Maximum number of records to insert in a single transaction
        
    Returns:
        Number of records successfully inserted
    """
    if not property_records:
        logger.info("No property records to insert")
        return 0
        
    # Group by property type for more efficient inserts
    grouped_properties = {}
    for record in property_records:
        prop_type = record.get('property_type_id')
        if not prop_type:
            logger.warning(f"Skipping record without property_type_id: {record}")
            continue
            
        if prop_type not in grouped_properties:
            grouped_properties[prop_type] = []
        grouped_properties[prop_type].append(record)
    
    # Insert each property type group in a single transaction
    total_inserted = 0
    total_groups = len(grouped_properties)
    
    logger.info(f"Inserting {len(property_records)} property records grouped by {total_groups} property types")
    
    for i, (prop_type, records) in enumerate(grouped_properties.items()):
        group_start_time = time.time()
        group_inserted = 0
        
        # Process in reasonable batch sizes
        for j in range(0, len(records), batch_size):
            batch = records[j:j+batch_size]
            batch_start_time = time.time()
            
            try:
                with safe_transaction() as conn:
                    # Build multi-value insert for this batch
                    columns = ['molecule_id', 'property_type_id', 'numeric_value',
                              'text_value', 'boolean_value', 'created_by']
                    values_list = []
                    params = []
                    
                    # Enhanced logging for debugging
                    logger.info(f"Preparing to insert {len(batch)} property records")
                    
                    for i, record in enumerate(batch):
                        # Validate required fields
                        if not record.get('molecule_id'):
                            logger.error(f"Record {i} missing molecule_id: {record}")
                            continue
                            
                        if not record.get('property_type_id'):
                            logger.error(f"Record {i} missing property_type_id: {record}")
                            continue
                        
                        values_list.append(f"(%s, %s, %s, %s, %s, %s)")
                        
                        # Log the actual values being inserted
                        molecule_id = record.get('molecule_id')
                        property_type_id = record.get('property_type_id')
                        numeric_value = record.get('numeric_value')
                        text_value = record.get('text_value')
                        boolean_value = record.get('boolean_value')
                        created_by = record.get('created_by')
                        
                        logger.debug(f"Record {i}: molecule_id={molecule_id}, property_type_id={property_type_id}, " +
                                    f"numeric_value={numeric_value}, text_value={text_value}, " +
                                    f"boolean_value={boolean_value}, created_by={created_by}")
                        
                        params.extend([
                            molecule_id,
                            property_type_id,
                            numeric_value,
                            text_value,
                            boolean_value,
                            created_by
                        ])
                    
                    if not values_list:
                        logger.warning("No valid property records to insert")
                        return 0
                    
                    # Execute the multi-value insert
                    query = f"""
                    INSERT INTO molecular_properties
                    ({', '.join(columns)})
                    VALUES {', '.join(values_list)}
                    RETURNING id, molecule_id, property_type_id
                    """
                    logger.info(f"Executing query: {query}")
                    logger.info(f"With params: {params}")
                    result = execute_query(query, params)
                    
                    # Verify the insertion by logging the returned IDs
                    if result:
                        logger.info(f"Successfully inserted {len(result)} property records")
                        for row in result:
                            logger.info(f"Inserted property: id={row.get('id')}, molecule_id={row.get('molecule_id')}, property_type_id={row.get('property_type_id')}")
                    else:
                        logger.warning("No property records were returned after insertion")
                    
                    batch_inserted = len(batch)
                    group_inserted += batch_inserted
                    total_inserted += batch_inserted
                    
                    batch_time = time.time() - batch_start_time
                    logger.debug(f"Inserted batch of {batch_inserted} records for property type {prop_type} "
                                f"({j+1}-{min(j+batch_size, len(records))}/{len(records)}) in {batch_time:.2f}s")
            
            except Exception as e:
                logger.error(f"Error inserting batch for property type {prop_type}: {str(e)}")
                # Log more details about the error
                import traceback
                logger.error(f"Error details: {traceback.format_exc()}")
                
                # Try to diagnose database-specific errors
                if "duplicate key value violates unique constraint" in str(e):
                    logger.error("Duplicate key violation - property may already exist")
                elif "foreign key constraint" in str(e):
                    logger.error("Foreign key constraint violation - property_type_id may not exist")
                elif "column" in str(e) and "does not exist" in str(e):
                    logger.error("Column does not exist - schema mismatch")
                
                # Continue with next batch despite errors
        
        group_time = time.time() - group_start_time
        logger.info(f"Inserted {group_inserted}/{len(records)} records for property type {prop_type} "
                   f"({i+1}/{total_groups}) in {group_time:.2f}s")
    
    logger.info(f"Total inserted: {total_inserted}/{len(property_records)} property records")
    return total_inserted


def resumable_batch_import(items: List[Any], process_func: callable, checkpoint_file: str, 
                          batch_size: int = 100) -> int:
    """
    Process items in batches with checkpoint-based resumability.
    
    This function implements the resumable batch import algorithm specified in
    DATABASE_POPULATION_ISSUES.md (Section 3.3.2). It processes items in batches
    and saves checkpoint information after each batch to allow resuming from the
    last successful batch in case of interruption.
    
    Args:
        items: List of items to process
        process_func: Function to call for each batch
        checkpoint_file: Path to the checkpoint file
        batch_size: Number of items to process in each batch
        
    Returns:
        Number of items successfully processed
    """
    if not items:
        logger.info("No items to process")
        return 0
        
    # Load checkpoint if exists
    checkpoint = {}
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
                logger.info(f"Loaded checkpoint from {checkpoint_file}: "
                           f"position={checkpoint.get('position', 0)}, "
                           f"processed={checkpoint.get('processed', 0)}")
        except Exception as e:
            logger.warning(f"Could not load checkpoint from {checkpoint_file}: {str(e)}")
    
    # Get starting position
    position = checkpoint.get('position', 0)
    processed_count = checkpoint.get('processed', 0)
    
    # Skip already processed items
    if position > 0:
        logger.info(f"Resuming from position {position} ({processed_count} items already processed)")
    
    # Process remaining items
    total_batches = math.ceil((len(items) - position) / batch_size)
    if total_batches <= 0:
        logger.info("All items have already been processed")
        return processed_count
        
    start_time = time.time()
    
    for batch_num in range(total_batches):
        batch_start = position + (batch_num * batch_size)
        batch_end = min(batch_start + batch_size, len(items))
        batch = items[batch_start:batch_end]
        
        batch_start_time = time.time()
        
        # Process this batch
        logger.info(f"Processing batch {batch_num+1}/{total_batches}, items {batch_start}-{batch_end}")
        try:
            batch_result = process_func(batch)
            batch_success = True
        except Exception as e:
            logger.error(f"Error processing batch {batch_num+1}: {str(e)}")
            batch_success = False
        
        # Update checkpoint only if batch was successful
        if batch_success:
            processed_count += len(batch)
            checkpoint = {
                'position': batch_end,
                'processed': processed_count,
                'last_updated': datetime.now().isoformat(),
                'last_batch': batch_num
            }
            
            # Save checkpoint
            try:
                with open(checkpoint_file, 'w') as f:
                    json.dump(checkpoint, f)
            except Exception as e:
                logger.error(f"Error saving checkpoint to {checkpoint_file}: {str(e)}")
        else:
            # If batch failed, don't update position so we'll retry this batch
            logger.warning(f"Batch {batch_num+1} failed, will retry on next run")
            break
                
        # Report progress
        batch_time = time.time() - batch_start_time
        elapsed = time.time() - start_time
        items_per_sec = processed_count / elapsed if elapsed > 0 else 0
        progress_pct = processed_count / len(items) * 100
        
        logger.info(f"Batch {batch_num+1} completed in {batch_time:.2f}s")
        logger.info(f"Progress: {processed_count}/{len(items)} items "
                   f"({progress_pct:.1f}%) at {items_per_sec:.1f} items/sec")
        
        # Estimate remaining time
        if batch_num < total_batches - 1:
            remaining_items = len(items) - batch_end
            est_remaining_time = remaining_items / items_per_sec if items_per_sec > 0 else 0
            est_completion = datetime.now() + datetime.timedelta(seconds=est_remaining_time)
            
            logger.info(f"Estimated remaining time: {est_remaining_time:.1f}s "
                       f"(completion around {est_completion.strftime('%H:%M:%S')})")
    
    total_time = time.time() - start_time
    logger.info(f"Batch processing completed: {processed_count}/{len(items)} items "
               f"processed in {total_time:.2f}s ({processed_count/total_time:.1f} items/sec)")
    
    return processed_count


def bulk_upsert_properties(property_records: List[Dict[str, Any]], batch_size: int = 1000) -> int:
    """
    Upsert multiple property records in optimized batches.
    
    This function is similar to bulk_insert_properties but uses an upsert operation
    (INSERT ... ON CONFLICT ... DO UPDATE) to update existing properties or insert
    new ones. This is useful for updating property values without having to check
    if they already exist.
    
    Args:
        property_records: List of property records to upsert
        batch_size: Maximum number of records to upsert in a single transaction
        
    Returns:
        Number of records successfully upserted
    """
    if not property_records:
        logger.info("No property records to upsert")
        return 0
        
    # Group by property type for more efficient upserts
    grouped_properties = {}
    for record in property_records:
        prop_type = record.get('property_type_id')
        if not prop_type:
            logger.warning(f"Skipping record without property_type_id: {record}")
            continue
            
        if prop_type not in grouped_properties:
            grouped_properties[prop_type] = []
        grouped_properties[prop_type].append(record)
    
    # Upsert each property type group in a single transaction
    total_upserted = 0
    total_groups = len(grouped_properties)
    
    logger.info(f"Upserting {len(property_records)} property records grouped by {total_groups} property types")
    
    for i, (prop_type, records) in enumerate(grouped_properties.items()):
        group_start_time = time.time()
        group_upserted = 0
        
        # Process in reasonable batch sizes
        for j in range(0, len(records), batch_size):
            batch = records[j:j+batch_size]
            batch_start_time = time.time()
            
            try:
                with safe_transaction() as conn:
                    # Build multi-value upsert for this batch
                    columns = ['molecule_id', 'property_type_id', 'numeric_value', 
                              'text_value', 'boolean_value', 'created_by']
                    values_list = []
                    params = []
                    
                    for record in batch:
                        values_list.append(f"(%s, %s, %s, %s, %s, %s)")
                        params.extend([
                            record.get('molecule_id'),
                            record.get('property_type_id'),
                            record.get('numeric_value'),
                            record.get('text_value'),
                            record.get('boolean_value'),
                            record.get('created_by')
                        ])
                    
                    # Execute the multi-value upsert
                    query = f"""
                    INSERT INTO molecular_properties 
                    ({', '.join(columns)})
                    VALUES {', '.join(values_list)}
                    ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                        numeric_value = EXCLUDED.numeric_value,
                        text_value = EXCLUDED.text_value,
                        boolean_value = EXCLUDED.boolean_value,
                        updated_at = NOW()
                    """
                    execute_query(query, params)
                    
                    batch_upserted = len(batch)
                    group_upserted += batch_upserted
                    total_upserted += batch_upserted
                    
                    batch_time = time.time() - batch_start_time
                    logger.debug(f"Upserted batch of {batch_upserted} records for property type {prop_type} "
                                f"({j+1}-{min(j+batch_size, len(records))}/{len(records)}) in {batch_time:.2f}s")
            
            except Exception as e:
                logger.error(f"Error upserting batch for property type {prop_type}: {str(e)}")
                # Continue with next batch despite errors
        
        group_time = time.time() - group_start_time
        logger.info(f"Upserted {group_upserted}/{len(records)} records for property type {prop_type} "
                   f"({i+1}/{total_groups}) in {group_time:.2f}s")
    
    logger.info(f"Total upserted: {total_upserted}/{len(property_records)} property records")
    return total_upserted


def batch_delete_properties(molecule_ids: List[UUID], property_type_ids: Optional[List[UUID]] = None,
                           batch_size: int = 1000) -> int:
    """
    Delete properties for multiple molecules in optimized batches.
    
    Args:
        molecule_ids: List of molecule UUIDs
        property_type_ids: Optional list of property type UUIDs (if None, delete all properties)
        batch_size: Maximum number of records to delete in a single transaction
        
    Returns:
        Number of property records deleted
    """
    if not molecule_ids:
        logger.info("No molecule IDs provided for property deletion")
        return 0
        
    total_deleted = 0
    total_batches = math.ceil(len(molecule_ids) / batch_size)
    
    logger.info(f"Deleting properties for {len(molecule_ids)} molecules "
               f"{'(all property types)' if property_type_ids is None else f'({len(property_type_ids)} property types)'}")
    
    for i in range(0, len(molecule_ids), batch_size):
        batch = molecule_ids[i:i+batch_size]
        batch_start_time = time.time()
        
        try:
            with safe_transaction() as conn:
                # Build the delete query
                molecule_placeholders = ', '.join(['%s'] * len(batch))
                
                if property_type_ids:
                    # Delete specific property types
                    property_placeholders = ', '.join(['%s'] * len(property_type_ids))
                    query = f"""
                    DELETE FROM molecular_properties
                    WHERE molecule_id IN ({molecule_placeholders})
                    AND property_type_id IN ({property_placeholders})
                    RETURNING molecule_id
                    """
                    params = batch + property_type_ids
                else:
                    # Delete all property types
                    query = f"""
                    DELETE FROM molecular_properties
                    WHERE molecule_id IN ({molecule_placeholders})
                    RETURNING molecule_id
                    """
                    params = batch
                
                # Execute the delete query
                result = execute_query(query, params)
                batch_deleted = len(result) if result else 0
                total_deleted += batch_deleted
                
                batch_time = time.time() - batch_start_time
                logger.info(f"Deleted {batch_deleted} property records for batch {i//batch_size + 1}/{total_batches} "
                           f"in {batch_time:.2f}s")
        
        except Exception as e:
            logger.error(f"Error deleting properties for batch {i//batch_size + 1}: {str(e)}")
            # Continue with next batch despite errors
    
    logger.info(f"Total deleted: {total_deleted} property records")
    return total_deleted


def batch_get_property_statistics(property_type_ids: List[UUID]) -> Dict[UUID, Dict[str, Any]]:
    """
    Get statistics for multiple property types in a single query.
    
    Args:
        property_type_ids: List of property type UUIDs
        
    Returns:
        Dictionary mapping property type UUIDs to statistics dictionaries
    """
    if not property_type_ids:
        logger.info("No property type IDs provided for statistics")
        return {}
        
    try:
        # Build the query
        placeholders = ', '.join(['%s'] * len(property_type_ids))
        query = f"""
        SELECT 
            property_type_id,
            COUNT(*) as count,
            MIN(numeric_value) as min_value,
            MAX(numeric_value) as max_value,
            AVG(numeric_value) as avg_value,
            STDDEV(numeric_value) as stddev_value
        FROM molecular_properties
        WHERE property_type_id IN ({placeholders})
        GROUP BY property_type_id
        """
        
        # Execute the query
        result = execute_query(query, property_type_ids)
        
        # Process the results
        stats = {}
        if result:
            for row in result:
                prop_id = row['property_type_id']
                stats[prop_id] = {
                    'count': row['count'],
                    'min_value': row['min_value'],
                    'max_value': row['max_value'],
                    'avg_value': row['avg_value'],
                    'stddev_value': row['stddev_value']
                }
        
        return stats
        
    except Exception as e:
        logger.error(f"Error getting property statistics: {str(e)}")
        return {}