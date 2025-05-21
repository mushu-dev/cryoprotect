"""
Bidirectional sync between Supabase and Convex.

This module provides utilities to synchronize data between Supabase and Convex,
ensuring data consistency across both systems.
"""

import os
import time
import json
import logging
import uuid
from typing import Dict, Any, List, Optional, Tuple, Callable, Union
from datetime import datetime, timedelta
import threading

# Import adapters
from .supabase_adapter import SupabaseDirectAdapter
from .enhanced_convex_adapter import ConvexAdapter

logger = logging.getLogger(__name__)

class ConvexSyncRecord:
    """Represents a record for synchronization tracking."""
    
    def __init__(self, table: str, record_id: str, source: str, last_sync: float, data: Dict):
        """
        Initialize a sync record.
        
        Args:
            table: Table name
            record_id: Record ID
            source: Source of the record ('supabase' or 'convex')
            last_sync: Last sync timestamp
            data: Record data
        """
        self.table = table
        self.record_id = record_id
        self.source = source
        self.last_sync = last_sync
        self.data = data
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'table': self.table,
            'record_id': self.record_id,
            'source': self.source,
            'last_sync': self.last_sync,
            'data': self.data
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'ConvexSyncRecord':
        """Create from dictionary."""
        return cls(
            table=data['table'],
            record_id=data['record_id'],
            source=data['source'],
            last_sync=data['last_sync'],
            data=data['data']
        )

class ConvexSyncManager:
    """
    Manager for bidirectional sync between Supabase and Convex.
    
    This class handles data synchronization between Supabase and Convex,
    ensuring consistency across both systems.
    """
    
    def __init__(
        self, 
        supabase_adapter: SupabaseDirectAdapter, 
        convex_adapter: ConvexAdapter,
        sync_interval: int = 60,
        sync_tables: Optional[List[str]] = None,
        conflict_resolution: str = 'last_modified'
    ):
        """
        Initialize the sync manager.
        
        Args:
            supabase_adapter: Supabase adapter
            convex_adapter: Convex adapter
            sync_interval: Sync interval in seconds
            sync_tables: List of tables to sync (None for all)
            conflict_resolution: Conflict resolution strategy
                ('last_modified', 'supabase_wins', 'convex_wins')
        """
        self.supabase = supabase_adapter
        self.convex = convex_adapter
        self.sync_interval = sync_interval
        self.sync_tables = sync_tables
        self.conflict_resolution = conflict_resolution
        
        # Sync state
        self.sync_state = {}
        self.last_sync_time = 0
        self.sync_in_progress = False
        self.sync_thread = None
        self.stop_sync = False
        
        # Callbacks
        self.on_sync_start = None
        self.on_sync_complete = None
        self.on_conflict = None
        self.on_sync_error = None
        
        # Modified since last sync
        self.modified_records = {
            'supabase': {},
            'convex': {}
        }
        
        # Initialize sync state
        self._init_sync_state()
    
    def _init_sync_state(self):
        """Initialize sync state tracking."""
        logger.info("Initializing sync state")
        
        # Get list of tables to sync
        tables = self.sync_tables
        
        if not tables:
            # Get tables from Supabase
            try:
                result = self.supabase.execute_query(
                    "SELECT table_name FROM information_schema.tables WHERE table_schema = 'public'"
                )
                tables = [row['table_name'] for row in result]
            except Exception as e:
                logger.error(f"Error getting tables from Supabase: {str(e)}")
                tables = []
        
        # Initialize sync state for each table
        for table in tables:
            self.sync_state[table] = {
                'last_sync': 0,
                'records': {}
            }
    
    def set_callbacks(
        self,
        on_sync_start: Optional[Callable] = None,
        on_sync_complete: Optional[Callable] = None,
        on_conflict: Optional[Callable[[str, str, Dict, Dict], Dict]] = None,
        on_sync_error: Optional[Callable[[Exception], None]] = None
    ):
        """
        Set sync callbacks.
        
        Args:
            on_sync_start: Called when sync starts
            on_sync_complete: Called when sync completes
            on_conflict: Called when a conflict is detected
                Function signature: (table, record_id, supabase_data, convex_data) -> resolved_data
            on_sync_error: Called when a sync error occurs
        """
        self.on_sync_start = on_sync_start
        self.on_sync_complete = on_sync_complete
        self.on_conflict = on_conflict
        self.on_sync_error = on_sync_error
    
    def start_sync(self, background: bool = True) -> bool:
        """
        Start sync process.
        
        Args:
            background: Run in background thread if True
            
        Returns:
            bool: True if sync started successfully
        """
        if self.sync_in_progress:
            logger.warning("Sync already in progress")
            return False
        
        if background:
            # Start sync in background thread
            self.stop_sync = False
            self.sync_thread = threading.Thread(target=self._sync_worker)
            self.sync_thread.daemon = True
            self.sync_thread.start()
            return True
        else:
            # Run sync immediately
            return self.sync_all()
    
    def stop_sync_thread(self):
        """Stop background sync thread."""
        self.stop_sync = True
        if self.sync_thread and self.sync_thread.is_alive():
            self.sync_thread.join(timeout=10)
            self.sync_thread = None
    
    def _sync_worker(self):
        """Background sync worker thread."""
        logger.info("Starting sync worker thread")
        
        while not self.stop_sync:
            try:
                # Run sync
                self.sync_all()
                
                # Sleep until next sync
                for _ in range(self.sync_interval):
                    if self.stop_sync:
                        break
                    time.sleep(1)
                    
            except Exception as e:
                logger.error(f"Error in sync worker thread: {str(e)}")
                if self.on_sync_error:
                    try:
                        self.on_sync_error(e)
                    except Exception as callback_error:
                        logger.error(f"Error in sync error callback: {str(callback_error)}")
                
                # Sleep after error to avoid spinning
                time.sleep(10)
        
        logger.info("Sync worker thread stopped")
    
    def sync_all(self) -> bool:
        """
        Synchronize all tables.
        
        Returns:
            bool: True if sync completed successfully
        """
        if self.sync_in_progress:
            logger.warning("Sync already in progress")
            return False
        
        self.sync_in_progress = True
        current_time = time.time()
        
        try:
            # Call sync start callback
            if self.on_sync_start:
                try:
                    self.on_sync_start()
                except Exception as e:
                    logger.error(f"Error in sync start callback: {str(e)}")
            
            logger.info(f"Starting sync at {datetime.now().isoformat()}")
            
            # Get tables to sync
            tables = list(self.sync_state.keys())
            
            # Sync each table
            success = True
            for table in tables:
                try:
                    if not self._sync_table(table):
                        success = False
                except Exception as e:
                    logger.error(f"Error syncing table {table}: {str(e)}")
                    if self.on_sync_error:
                        try:
                            self.on_sync_error(e)
                        except Exception as callback_error:
                            logger.error(f"Error in sync error callback: {str(callback_error)}")
                    success = False
            
            # Update last sync time
            self.last_sync_time = current_time
            
            # Call sync complete callback
            if self.on_sync_complete:
                try:
                    self.on_sync_complete()
                except Exception as e:
                    logger.error(f"Error in sync complete callback: {str(e)}")
            
            logger.info(f"Sync completed at {datetime.now().isoformat()}")
            
            return success
            
        finally:
            self.sync_in_progress = False
    
    def _sync_table(self, table: str) -> bool:
        """
        Synchronize a single table.
        
        Args:
            table: Table name
            
        Returns:
            bool: True if sync completed successfully
        """
        logger.info(f"Syncing table: {table}")
        
        try:
            # Get last sync time for table
            last_sync = self.sync_state[table]['last_sync']
            
            # Get records from Supabase modified since last sync
            supabase_records = self._get_modified_records_from_supabase(table, last_sync)
            
            # Get records from Convex modified since last sync
            convex_records = self._get_modified_records_from_convex(table, last_sync)
            
            # Detect conflicts and resolve
            conflicts = self._detect_conflicts(table, supabase_records, convex_records)
            
            # Apply non-conflicting changes
            self._apply_non_conflicting_changes(table, supabase_records, convex_records, conflicts)
            
            # Resolve and apply conflicting changes
            self._resolve_conflicts(table, conflicts)
            
            # Update sync state
            self.sync_state[table]['last_sync'] = time.time()
            
            return True
            
        except Exception as e:
            logger.error(f"Error syncing table {table}: {str(e)}")
            return False
    
    def _get_modified_records_from_supabase(self, table: str, since: float) -> Dict[str, Dict]:
        """
        Get records modified in Supabase since a given time.
        
        Args:
            table: Table name
            since: Timestamp
            
        Returns:
            Dict mapping record IDs to record data
        """
        # Convert timestamp to ISO format
        since_date = datetime.fromtimestamp(since).isoformat()
        
        try:
            # Get records modified since last sync
            # This assumes tables have 'updated_at' column
            result = self.supabase.table(table) \
                .select('*') \
                .gte('updated_at', since_date) \
                .execute()
            
            records = {}
            for record in result.data:
                record_id = str(record.get('id'))
                if record_id:
                    records[record_id] = record
            
            logger.info(f"Found {len(records)} modified records in Supabase table {table}")
            return records
            
        except Exception as e:
            logger.error(f"Error getting modified records from Supabase table {table}: {str(e)}")
            return {}
    
    def _get_modified_records_from_convex(self, table: str, since: float) -> Dict[str, Dict]:
        """
        Get records modified in Convex since a given time.
        
        Args:
            table: Table name
            since: Timestamp
            
        Returns:
            Dict mapping record IDs to record data
        """
        # Convert timestamp to ISO format
        since_date = datetime.fromtimestamp(since).isoformat()
        
        try:
            # Get records modified since last sync
            # This assumes tables have 'updatedAt' field
            result = self.convex.execute_query(f"api.sync.getModifiedRecords", {
                'tableName': table,
                'since': since_date
            })
            
            records = {}
            if result and isinstance(result, list):
                for record in result:
                    record_id = str(record.get('id'))
                    if record_id:
                        records[record_id] = record
            
            logger.info(f"Found {len(records)} modified records in Convex table {table}")
            return records
            
        except Exception as e:
            logger.error(f"Error getting modified records from Convex table {table}: {str(e)}")
            return {}
    
    def _detect_conflicts(
        self, 
        table: str,
        supabase_records: Dict[str, Dict], 
        convex_records: Dict[str, Dict]
    ) -> Dict[str, Tuple[Dict, Dict]]:
        """
        Detect conflicts between Supabase and Convex records.
        
        Args:
            table: Table name
            supabase_records: Records from Supabase
            convex_records: Records from Convex
            
        Returns:
            Dict mapping record IDs to tuples of (supabase_record, convex_record)
        """
        conflicts = {}
        
        # Find records modified in both systems
        for record_id in set(supabase_records.keys()) & set(convex_records.keys()):
            supabase_record = supabase_records[record_id]
            convex_record = convex_records[record_id]
            
            # Check if records are different
            if self._records_differ(supabase_record, convex_record):
                conflicts[record_id] = (supabase_record, convex_record)
        
        logger.info(f"Detected {len(conflicts)} conflicts in table {table}")
        return conflicts
    
    def _records_differ(self, record1: Dict, record2: Dict) -> bool:
        """
        Check if two records differ in their content.
        
        Args:
            record1: First record
            record2: Second record
            
        Returns:
            bool: True if records differ
        """
        # Fields to ignore in comparison
        ignore_fields = ['updated_at', 'updatedAt', '_creationTime', '_id']
        
        for key in set(record1.keys()) | set(record2.keys()):
            if key in ignore_fields:
                continue
                
            if key not in record1 or key not in record2:
                return True
                
            if record1[key] != record2[key]:
                return True
        
        return False
    
    def _apply_non_conflicting_changes(
        self,
        table: str,
        supabase_records: Dict[str, Dict],
        convex_records: Dict[str, Dict],
        conflicts: Dict[str, Tuple[Dict, Dict]]
    ) -> None:
        """
        Apply non-conflicting changes.
        
        Args:
            table: Table name
            supabase_records: Records from Supabase
            convex_records: Records from Convex
            conflicts: Detected conflicts
        """
        # Apply Supabase changes to Convex
        for record_id, record in supabase_records.items():
            if record_id not in conflicts and record_id not in convex_records:
                # Record exists in Supabase but not in Convex
                self._create_convex_record(table, record)
                
        # Apply Convex changes to Supabase
        for record_id, record in convex_records.items():
            if record_id not in conflicts and record_id not in supabase_records:
                # Record exists in Convex but not in Supabase
                self._create_supabase_record(table, record)
    
    def _resolve_conflicts(self, table: str, conflicts: Dict[str, Tuple[Dict, Dict]]) -> None:
        """
        Resolve conflicts.
        
        Args:
            table: Table name
            conflicts: Detected conflicts
        """
        for record_id, (supabase_record, convex_record) in conflicts.items():
            resolved_record = None
            
            # Use custom conflict resolution if provided
            if self.on_conflict:
                try:
                    resolved_record = self.on_conflict(table, record_id, supabase_record, convex_record)
                except Exception as e:
                    logger.error(f"Error in conflict resolution callback: {str(e)}")
            
            if not resolved_record:
                # Use configured conflict resolution strategy
                if self.conflict_resolution == 'last_modified':
                    # Use record with later modification time
                    supabase_time = self._get_record_modified_time(supabase_record)
                    convex_time = self._get_record_modified_time(convex_record)
                    
                    if supabase_time >= convex_time:
                        resolved_record = supabase_record
                    else:
                        resolved_record = convex_record
                        
                elif self.conflict_resolution == 'supabase_wins':
                    resolved_record = supabase_record
                else:  # convex_wins
                    resolved_record = convex_record
            
            # Apply resolved record to both systems
            if resolved_record:
                self._update_supabase_record(table, resolved_record)
                self._update_convex_record(table, resolved_record)
    
    def _get_record_modified_time(self, record: Dict) -> float:
        """
        Get modification time from record.
        
        Args:
            record: Record data
            
        Returns:
            float: Modification timestamp
        """
        # Try different field names for modification time
        for field in ['updated_at', 'updatedAt', '_updatedAt', 'modified_at', 'modifiedAt']:
            if field in record:
                value = record[field]
                
                # Handle different formats
                if isinstance(value, (int, float)):
                    return float(value)
                elif isinstance(value, str):
                    try:
                        dt = datetime.fromisoformat(value.replace('Z', '+00:00'))
                        return dt.timestamp()
                    except (ValueError, TypeError):
                        pass
        
        # Default to current time if no field found
        return time.time()
    
    def _create_convex_record(self, table: str, record: Dict) -> None:
        """
        Create a record in Convex.
        
        Args:
            table: Table name
            record: Record data
        """
        try:
            # Prepare record for Convex
            convex_record = self._prepare_record_for_convex(record)
            
            # Create record in Convex
            self.convex.execute_query(f"api.sync.createRecord", {
                'tableName': table,
                'record': convex_record
            })
            
            logger.info(f"Created record in Convex table {table}: {record.get('id')}")
            
        except Exception as e:
            logger.error(f"Error creating record in Convex table {table}: {str(e)}")
    
    def _create_supabase_record(self, table: str, record: Dict) -> None:
        """
        Create a record in Supabase.
        
        Args:
            table: Table name
            record: Record data
        """
        try:
            # Prepare record for Supabase
            supabase_record = self._prepare_record_for_supabase(record)
            
            # Create record in Supabase
            self.supabase.table(table).insert(supabase_record).execute()
            
            logger.info(f"Created record in Supabase table {table}: {record.get('id')}")
            
        except Exception as e:
            logger.error(f"Error creating record in Supabase table {table}: {str(e)}")
    
    def _update_convex_record(self, table: str, record: Dict) -> None:
        """
        Update a record in Convex.
        
        Args:
            table: Table name
            record: Record data
        """
        try:
            # Get record ID
            record_id = record.get('id')
            if not record_id:
                logger.error(f"Cannot update record in Convex: Missing ID")
                return
            
            # Prepare record for Convex
            convex_record = self._prepare_record_for_convex(record)
            
            # Update record in Convex
            self.convex.execute_query(f"api.sync.updateRecord", {
                'tableName': table,
                'recordId': record_id,
                'record': convex_record
            })
            
            logger.info(f"Updated record in Convex table {table}: {record_id}")
            
        except Exception as e:
            logger.error(f"Error updating record in Convex table {table}: {str(e)}")
    
    def _update_supabase_record(self, table: str, record: Dict) -> None:
        """
        Update a record in Supabase.
        
        Args:
            table: Table name
            record: Record data
        """
        try:
            # Get record ID
            record_id = record.get('id')
            if not record_id:
                logger.error(f"Cannot update record in Supabase: Missing ID")
                return
            
            # Prepare record for Supabase
            supabase_record = self._prepare_record_for_supabase(record)
            
            # Update record in Supabase
            self.supabase.table(table) \
                .update(supabase_record) \
                .eq('id', record_id) \
                .execute()
            
            logger.info(f"Updated record in Supabase table {table}: {record_id}")
            
        except Exception as e:
            logger.error(f"Error updating record in Supabase table {table}: {str(e)}")
    
    def _prepare_record_for_convex(self, record: Dict) -> Dict:
        """
        Prepare record for Convex.
        
        Args:
            record: Record data
            
        Returns:
            Dict: Record prepared for Convex
        """
        # Create a copy to avoid modifying the original
        result = record.copy()
        
        # Convert field names to camelCase if needed
        field_mapping = {
            'updated_at': 'updatedAt',
            'created_at': 'createdAt'
        }
        
        for snake_case, camel_case in field_mapping.items():
            if snake_case in result and camel_case not in result:
                result[camel_case] = result[snake_case]
        
        # Remove Convex-specific fields
        convex_fields = ['_id', '_creationTime']
        for field in convex_fields:
            if field in result:
                del result[field]
        
        return result
    
    def _prepare_record_for_supabase(self, record: Dict) -> Dict:
        """
        Prepare record for Supabase.
        
        Args:
            record: Record data
            
        Returns:
            Dict: Record prepared for Supabase
        """
        # Create a copy to avoid modifying the original
        result = record.copy()
        
        # Convert field names to snake_case if needed
        field_mapping = {
            'updatedAt': 'updated_at',
            'createdAt': 'created_at'
        }
        
        for camel_case, snake_case in field_mapping.items():
            if camel_case in result and snake_case not in result:
                result[snake_case] = result[camel_case]
        
        # Remove Convex-specific fields
        convex_fields = ['_id', '_creationTime']
        for field in convex_fields:
            if field in result:
                del result[field]
        
        return result
    
    def force_sync_record(self, table: str, record_id: str, target: str = 'both') -> bool:
        """
        Force sync a specific record.
        
        Args:
            table: Table name
            record_id: Record ID
            target: Sync target ('both', 'supabase', 'convex')
            
        Returns:
            bool: True if sync successful
        """
        try:
            # Get record from Supabase
            supabase_record = None
            if target in ('both', 'convex'):
                result = self.supabase.table(table) \
                    .select('*') \
                    .eq('id', record_id) \
                    .execute()
                
                if result.data and len(result.data) > 0:
                    supabase_record = result.data[0]
            
            # Get record from Convex
            convex_record = None
            if target in ('both', 'supabase'):
                result = self.convex.execute_query(f"api.sync.getRecord", {
                    'tableName': table,
                    'recordId': record_id
                })
                
                if result:
                    convex_record = result
            
            # Sync based on target
            if target == 'both':
                # Check if both records exist
                if supabase_record and convex_record:
                    # Detect and resolve conflict
                    if self._records_differ(supabase_record, convex_record):
                        conflicts = {record_id: (supabase_record, convex_record)}
                        self._resolve_conflicts(table, conflicts)
                elif supabase_record:
                    # Create in Convex
                    self._create_convex_record(table, supabase_record)
                elif convex_record:
                    # Create in Supabase
                    self._create_supabase_record(table, convex_record)
                else:
                    logger.warning(f"Record {record_id} not found in either system")
                    return False
                    
            elif target == 'convex' and supabase_record:
                # Update Convex from Supabase
                if convex_record:
                    self._update_convex_record(table, supabase_record)
                else:
                    self._create_convex_record(table, supabase_record)
                    
            elif target == 'supabase' and convex_record:
                # Update Supabase from Convex
                if supabase_record:
                    self._update_supabase_record(table, convex_record)
                else:
                    self._create_supabase_record(table, convex_record)
                    
            else:
                logger.warning(f"Cannot sync record {record_id}: No source record found")
                return False
            
            return True
            
        except Exception as e:
            logger.error(f"Error forcing sync for record {record_id}: {str(e)}")
            if self.on_sync_error:
                try:
                    self.on_sync_error(e)
                except Exception as callback_error:
                    logger.error(f"Error in sync error callback: {str(callback_error)}")
            return False

def create_sync_manager(
    supabase_adapter: SupabaseDirectAdapter,
    convex_adapter: ConvexAdapter,
    config: Optional[Dict[str, Any]] = None
) -> ConvexSyncManager:
    """
    Create a ConvexSyncManager instance.
    
    Args:
        supabase_adapter: Supabase adapter
        convex_adapter: Convex adapter
        config: Optional configuration dict
        
    Returns:
        ConvexSyncManager: Configured ConvexSyncManager instance
    """
    # Default configuration
    default_config = {
        'sync_interval': int(os.environ.get('CONVEX_SYNC_INTERVAL', '60')),
        'sync_tables': os.environ.get('CONVEX_SYNC_TABLES', '').split(',') if os.environ.get('CONVEX_SYNC_TABLES') else None,
        'conflict_resolution': os.environ.get('CONVEX_CONFLICT_RESOLUTION', 'last_modified')
    }
    
    # Merge with provided config
    if config:
        default_config.update(config)
    
    # Create sync manager
    return ConvexSyncManager(
        supabase_adapter=supabase_adapter,
        convex_adapter=convex_adapter,
        sync_interval=default_config['sync_interval'],
        sync_tables=default_config['sync_tables'],
        conflict_resolution=default_config['conflict_resolution']
    )