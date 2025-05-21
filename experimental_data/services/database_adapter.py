#!/usr/bin/env python3
"""
Database Adapter for CryoProtect Enhanced Experimental Data System.

This module provides a database adapter for the experimental data system,
which handles database connection, query execution, and data transformation.
"""

from typing import Dict, List, Any, Optional, Union, Tuple, Callable
import uuid
from datetime import datetime
import json
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

class DatabaseAdapter:
    """Base class for database adapters."""
    
    async def create(self, table: str, data: Dict[str, Any]) -> str:
        """
        Create a new record.
        
        Args:
            table: Table name
            data: Record data
            
        Returns:
            ID of the created record
        """
        raise NotImplementedError("Subclasses must implement create()")
    
    async def get(self, table: str, id: str) -> Optional[Dict[str, Any]]:
        """
        Get a record by ID.
        
        Args:
            table: Table name
            id: Record ID
            
        Returns:
            Record data if found, None otherwise
        """
        raise NotImplementedError("Subclasses must implement get()")
    
    async def list(
        self, 
        table: str,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: int = 20,
        sort_by: str = 'created_at',
        sort_order: str = 'desc'
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List records with pagination and filtering.
        
        Args:
            table: Table name
            filters: Filters to apply
            offset: Offset for pagination
            limit: Limit for pagination
            sort_by: Field to sort by
            sort_order: Sort order ('asc' or 'desc')
            
        Returns:
            Tuple of (records, total_count)
        """
        raise NotImplementedError("Subclasses must implement list()")
    
    async def update(self, table: str, id: str, data: Dict[str, Any]) -> bool:
        """
        Update a record.
        
        Args:
            table: Table name
            id: Record ID
            data: Updated data
            
        Returns:
            True if updated, False otherwise
        """
        raise NotImplementedError("Subclasses must implement update()")
    
    async def delete(self, table: str, id: str) -> bool:
        """
        Delete a record.
        
        Args:
            table: Table name
            id: Record ID
            
        Returns:
            True if deleted, False otherwise
        """
        raise NotImplementedError("Subclasses must implement delete()")
    
    async def count(self, table: str, filters: Optional[Dict[str, Any]] = None) -> int:
        """
        Count records matching filters.
        
        Args:
            table: Table name
            filters: Filters to apply
            
        Returns:
            Count of matching records
        """
        raise NotImplementedError("Subclasses must implement count()")
    
    async def execute_query(self, query: str, params: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
        """
        Execute a custom query.
        
        Args:
            query: SQL query
            params: Query parameters
            
        Returns:
            Query results
        """
        raise NotImplementedError("Subclasses must implement execute_query()")
    
    async def begin_transaction(self) -> None:
        """Begin a transaction."""
        raise NotImplementedError("Subclasses must implement begin_transaction()")
    
    async def commit_transaction(self) -> None:
        """Commit the current transaction."""
        raise NotImplementedError("Subclasses must implement commit_transaction()")
    
    async def rollback_transaction(self) -> None:
        """Roll back the current transaction."""
        raise NotImplementedError("Subclasses must implement rollback_transaction()")

class SupabaseAdapter(DatabaseAdapter):
    """Database adapter for Supabase."""
    
    def __init__(self, supabase_client=None):
        """
        Initialize Supabase adapter.
        
        Args:
            supabase_client: Supabase client instance (optional)
        """
        self.supabase = supabase_client
        self._transaction_in_progress = False
    
    async def _ensure_client(self):
        """Ensure a Supabase client is available."""
        if not self.supabase:
            from service_role_helper import get_supabase_client
            self.supabase = get_supabase_client()
    
    async def create(self, table: str, data: Dict[str, Any]) -> str:
        """
        Create a new record.
        
        Args:
            table: Table name
            data: Record data
            
        Returns:
            ID of the created record
        """
        await self._ensure_client()
        
        # Ensure data has an ID
        if 'id' not in data:
            data['id'] = str(uuid.uuid4())
        
        # Ensure created_at and updated_at
        if 'created_at' not in data:
            data['created_at'] = datetime.now().isoformat()
        
        if 'updated_at' not in data:
            data['updated_at'] = datetime.now().isoformat()
        
        # Create record
        result = await self.supabase.table(table).insert(data).execute()
        
        if not hasattr(result, 'data') or not result.data:
            raise Exception(f"Failed to create record in {table}")
        
        return result.data[0]['id']
    
    async def get(self, table: str, id: str) -> Optional[Dict[str, Any]]:
        """
        Get a record by ID.
        
        Args:
            table: Table name
            id: Record ID
            
        Returns:
            Record data if found, None otherwise
        """
        await self._ensure_client()
        
        result = await self.supabase.table(table).select('*').eq('id', id).single().execute()
        
        if not hasattr(result, 'data') or not result.data:
            return None
        
        return result.data
    
    async def list(
        self, 
        table: str,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: int = 20,
        sort_by: str = 'created_at',
        sort_order: str = 'desc'
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List records with pagination and filtering.
        
        Args:
            table: Table name
            filters: Filters to apply
            offset: Offset for pagination
            limit: Limit for pagination
            sort_by: Field to sort by
            sort_order: Sort order ('asc' or 'desc')
            
        Returns:
            Tuple of (records, total_count)
        """
        await self._ensure_client()
        
        # Build query
        query = self.supabase.table(table).select('*', count='exact')
        
        # Apply filters
        if filters:
            for field, value in filters.items():
                # Handle special operator suffixes
                if '__' in field:
                    field_name, operator = field.split('__', 1)
                    
                    if operator == 'eq':
                        query = query.eq(field_name, value)
                    elif operator == 'neq':
                        query = query.neq(field_name, value)
                    elif operator == 'gt':
                        query = query.gt(field_name, value)
                    elif operator == 'gte':
                        query = query.gte(field_name, value)
                    elif operator == 'lt':
                        query = query.lt(field_name, value)
                    elif operator == 'lte':
                        query = query.lte(field_name, value)
                    elif operator == 'like':
                        query = query.like(field_name, f'%{value}%')
                    elif operator == 'ilike':
                        query = query.ilike(field_name, f'%{value}%')
                    elif operator == 'in':
                        query = query.in_(field_name, value)
                    elif operator == 'is':
                        if value is None:
                            query = query.is_(field_name, 'null')
                        else:
                            query = query.is_(field_name, 'not.null')
                else:
                    query = query.eq(field, value)
        
        # Apply sorting
        if sort_order.lower() == 'desc':
            query = query.order(sort_by, desc=True)
        else:
            query = query.order(sort_by)
        
        # Apply pagination
        query = query.range(offset, offset + limit - 1)
        
        # Execute query
        result = await query.execute()
        
        if not hasattr(result, 'data'):
            return [], 0
        
        count = result.count if hasattr(result, 'count') else len(result.data)
        
        return result.data, count
    
    async def update(self, table: str, id: str, data: Dict[str, Any]) -> bool:
        """
        Update a record.
        
        Args:
            table: Table name
            id: Record ID
            data: Updated data
            
        Returns:
            True if updated, False otherwise
        """
        await self._ensure_client()
        
        # Ensure updated_at
        if 'updated_at' not in data:
            data['updated_at'] = datetime.now().isoformat()
        
        # Remove id from data if present
        if 'id' in data:
            del data['id']
        
        # Update record
        result = await self.supabase.table(table).update(data).eq('id', id).execute()
        
        if not hasattr(result, 'data') or not result.data:
            return False
        
        return True
    
    async def delete(self, table: str, id: str) -> bool:
        """
        Delete a record.
        
        Args:
            table: Table name
            id: Record ID
            
        Returns:
            True if deleted, False otherwise
        """
        await self._ensure_client()
        
        result = await self.supabase.table(table).delete().eq('id', id).execute()
        
        if not hasattr(result, 'data') or not result.data:
            return False
        
        return True
    
    async def count(self, table: str, filters: Optional[Dict[str, Any]] = None) -> int:
        """
        Count records matching filters.
        
        Args:
            table: Table name
            filters: Filters to apply
            
        Returns:
            Count of matching records
        """
        await self._ensure_client()
        
        # Build query
        query = self.supabase.table(table).select('id', count='exact')
        
        # Apply filters
        if filters:
            for field, value in filters.items():
                # Handle special operator suffixes
                if '__' in field:
                    field_name, operator = field.split('__', 1)
                    
                    if operator == 'eq':
                        query = query.eq(field_name, value)
                    elif operator == 'neq':
                        query = query.neq(field_name, value)
                    elif operator == 'gt':
                        query = query.gt(field_name, value)
                    elif operator == 'gte':
                        query = query.gte(field_name, value)
                    elif operator == 'lt':
                        query = query.lt(field_name, value)
                    elif operator == 'lte':
                        query = query.lte(field_name, value)
                    elif operator == 'like':
                        query = query.like(field_name, f'%{value}%')
                    elif operator == 'ilike':
                        query = query.ilike(field_name, f'%{value}%')
                    elif operator == 'in':
                        query = query.in_(field_name, value)
                    elif operator == 'is':
                        if value is None:
                            query = query.is_(field_name, 'null')
                        else:
                            query = query.is_(field_name, 'not.null')
                else:
                    query = query.eq(field, value)
        
        # Execute query
        result = await query.execute()
        
        if not hasattr(result, 'count'):
            return 0
        
        return result.count
    
    async def execute_query(
        self, 
        query: str, 
        params: Optional[Dict[str, Any]] = None
    ) -> List[Dict[str, Any]]:
        """
        Execute a custom query.
        
        Args:
            query: SQL query
            params: Query parameters
            
        Returns:
            Query results
        """
        from supabase_direct import execute_sql_query
        
        result = await execute_sql_query(query, params)
        return result
    
    async def begin_transaction(self) -> None:
        """Begin a transaction."""
        if self._transaction_in_progress:
            raise Exception("Transaction already in progress")
        
        self._transaction_in_progress = True
        await self.execute_query("BEGIN")
    
    async def commit_transaction(self) -> None:
        """Commit the current transaction."""
        if not self._transaction_in_progress:
            raise Exception("No transaction in progress")
        
        await self.execute_query("COMMIT")
        self._transaction_in_progress = False
    
    async def rollback_transaction(self) -> None:
        """Roll back the current transaction."""
        if not self._transaction_in_progress:
            raise Exception("No transaction in progress")
        
        await self.execute_query("ROLLBACK")
        self._transaction_in_progress = False

class MockDatabaseAdapter(DatabaseAdapter):
    """Mock database adapter for testing."""
    
    def __init__(self):
        """Initialize mock database."""
        self.data = {
            'experiments': {},
            'experiment_results': {},
            'protocols': {},
            'tissue_types': {},
            'time_series': {},
            'time_series_data': {},
            'validation_rules': {}
        }
        self._transaction_in_progress = False
        self._transaction_backup = None
    
    async def create(self, table: str, data: Dict[str, Any]) -> str:
        """
        Create a new record.
        
        Args:
            table: Table name
            data: Record data
            
        Returns:
            ID of the created record
        """
        # Ensure table exists
        if table not in self.data:
            self.data[table] = {}
        
        # Ensure data has an ID
        if 'id' not in data:
            data['id'] = str(uuid.uuid4())
        
        # Ensure created_at and updated_at
        if 'created_at' not in data:
            data['created_at'] = datetime.now().isoformat()
        
        if 'updated_at' not in data:
            data['updated_at'] = datetime.now().isoformat()
        
        # Store record
        self.data[table][data['id']] = data.copy()
        
        return data['id']
    
    async def get(self, table: str, id: str) -> Optional[Dict[str, Any]]:
        """
        Get a record by ID.
        
        Args:
            table: Table name
            id: Record ID
            
        Returns:
            Record data if found, None otherwise
        """
        if table not in self.data or id not in self.data[table]:
            return None
        
        return self.data[table][id].copy()
    
    async def list(
        self, 
        table: str,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: int = 20,
        sort_by: str = 'created_at',
        sort_order: str = 'desc'
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List records with pagination and filtering.
        
        Args:
            table: Table name
            filters: Filters to apply
            offset: Offset for pagination
            limit: Limit for pagination
            sort_by: Field to sort by
            sort_order: Sort order ('asc' or 'desc')
            
        Returns:
            Tuple of (records, total_count)
        """
        if table not in self.data:
            return [], 0
        
        # Get all records
        records = list(self.data[table].values())
        
        # Apply filters
        if filters:
            filtered_records = []
            for record in records:
                match = True
                
                for field, value in filters.items():
                    # Handle special operator suffixes
                    if '__' in field:
                        field_name, operator = field.split('__', 1)
                        
                        if field_name not in record:
                            match = False
                            break
                        
                        record_value = record[field_name]
                        
                        if operator == 'eq' and record_value != value:
                            match = False
                            break
                        elif operator == 'neq' and record_value == value:
                            match = False
                            break
                        elif operator == 'gt' and record_value <= value:
                            match = False
                            break
                        elif operator == 'gte' and record_value < value:
                            match = False
                            break
                        elif operator == 'lt' and record_value >= value:
                            match = False
                            break
                        elif operator == 'lte' and record_value > value:
                            match = False
                            break
                        elif operator == 'like':
                            if not isinstance(record_value, str) or value not in record_value:
                                match = False
                                break
                        elif operator == 'ilike':
                            if not isinstance(record_value, str) or value.lower() not in record_value.lower():
                                match = False
                                break
                        elif operator == 'in' and record_value not in value:
                            match = False
                            break
                        elif operator == 'is':
                            if value is None and record_value is not None:
                                match = False
                                break
                            elif value is not None and record_value is None:
                                match = False
                                break
                    else:
                        if field not in record or record[field] != value:
                            match = False
                            break
                
                if match:
                    filtered_records.append(record)
            
            records = filtered_records
        
        # Sort records
        reverse = sort_order.lower() == 'desc'
        records.sort(key=lambda r: r.get(sort_by) if sort_by in r else None, reverse=reverse)
        
        # Apply pagination
        paginated_records = records[offset:offset + limit]
        
        return [r.copy() for r in paginated_records], len(records)
    
    async def update(self, table: str, id: str, data: Dict[str, Any]) -> bool:
        """
        Update a record.
        
        Args:
            table: Table name
            id: Record ID
            data: Updated data
            
        Returns:
            True if updated, False otherwise
        """
        if table not in self.data or id not in self.data[table]:
            return False
        
        # Get existing record
        record = self.data[table][id]
        
        # Update fields
        for key, value in data.items():
            if key != 'id':  # Don't update ID
                record[key] = value
        
        # Update updated_at
        record['updated_at'] = datetime.now().isoformat()
        
        return True
    
    async def delete(self, table: str, id: str) -> bool:
        """
        Delete a record.
        
        Args:
            table: Table name
            id: Record ID
            
        Returns:
            True if deleted, False otherwise
        """
        if table not in self.data or id not in self.data[table]:
            return False
        
        del self.data[table][id]
        return True
    
    async def count(self, table: str, filters: Optional[Dict[str, Any]] = None) -> int:
        """
        Count records matching filters.
        
        Args:
            table: Table name
            filters: Filters to apply
            
        Returns:
            Count of matching records
        """
        if table not in self.data:
            return 0
        
        # Get all records
        records = list(self.data[table].values())
        
        # Apply filters
        if filters:
            filtered_records = []
            for record in records:
                match = True
                
                for field, value in filters.items():
                    # Handle special operator suffixes
                    if '__' in field:
                        field_name, operator = field.split('__', 1)
                        
                        if field_name not in record:
                            match = False
                            break
                        
                        record_value = record[field_name]
                        
                        if operator == 'eq' and record_value != value:
                            match = False
                            break
                        elif operator == 'neq' and record_value == value:
                            match = False
                            break
                        elif operator == 'gt' and record_value <= value:
                            match = False
                            break
                        elif operator == 'gte' and record_value < value:
                            match = False
                            break
                        elif operator == 'lt' and record_value >= value:
                            match = False
                            break
                        elif operator == 'lte' and record_value > value:
                            match = False
                            break
                        elif operator == 'like':
                            if not isinstance(record_value, str) or value not in record_value:
                                match = False
                                break
                        elif operator == 'ilike':
                            if not isinstance(record_value, str) or value.lower() not in record_value.lower():
                                match = False
                                break
                        elif operator == 'in' and record_value not in value:
                            match = False
                            break
                        elif operator == 'is':
                            if value is None and record_value is not None:
                                match = False
                                break
                            elif value is not None and record_value is None:
                                match = False
                                break
                    else:
                        if field not in record or record[field] != value:
                            match = False
                            break
                
                if match:
                    filtered_records.append(record)
            
            records = filtered_records
        
        return len(records)
    
    async def execute_query(
        self, 
        query: str, 
        params: Optional[Dict[str, Any]] = None
    ) -> List[Dict[str, Any]]:
        """
        Execute a custom query.
        
        Args:
            query: SQL query
            params: Query parameters
            
        Returns:
            Query results
        """
        # This is a mock adapter, so it doesn't support raw SQL queries
        raise NotImplementedError("Mock adapter does not support raw SQL queries")
    
    async def begin_transaction(self) -> None:
        """Begin a transaction."""
        if self._transaction_in_progress:
            raise Exception("Transaction already in progress")
        
        # Create a deep copy of the current data state
        import copy
        self._transaction_backup = copy.deepcopy(self.data)
        self._transaction_in_progress = True
    
    async def commit_transaction(self) -> None:
        """Commit the current transaction."""
        if not self._transaction_in_progress:
            raise Exception("No transaction in progress")
        
        # Clear the backup and end the transaction
        self._transaction_backup = None
        self._transaction_in_progress = False
    
    async def rollback_transaction(self) -> None:
        """Roll back the current transaction."""
        if not self._transaction_in_progress:
            raise Exception("No transaction in progress")
        
        # Restore the data from the backup
        self.data = self._transaction_backup
        self._transaction_backup = None
        self._transaction_in_progress = False

def create_database_adapter() -> DatabaseAdapter:
    """
    Create a database adapter based on environment configuration.
    
    Returns:
        Database adapter instance
    """
    adapter_type = os.getenv('DB_ADAPTER_TYPE', 'supabase')
    
    if adapter_type.lower() == 'mock':
        return MockDatabaseAdapter()
    else:
        return SupabaseAdapter()