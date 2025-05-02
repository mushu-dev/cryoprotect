"""
Supabase mock implementation.

This module provides a mock implementation of the Supabase client
for testing without connecting to an actual Supabase instance.
"""

import pytest
from unittest.mock import patch
from contextlib import contextmanager
from typing import Dict, List, Any, Optional, Union, Callable, Generator

class MockQueryResult:
    """Mock Supabase query result."""
    
    def __init__(self, data: List[Dict] = None, error: Any = None):
        self.data = data or []
        self.error = error
        
class MockQuery:
    """Mock Supabase query builder."""
    
    def __init__(self, table: str, data: List[Dict]):
        self.table = table
        self.data = data
        self.filters = []
        self.selected_columns = None
        self._order_by = None
        self._limit = None
        
    def select(self, columns: str = '*') -> 'MockQuery':
        """Mock select method."""
        self.selected_columns = columns
        return self
        
    def eq(self, column: str, value: Any) -> 'MockQuery':
        """Mock equality filter."""
        self.filters.append(('eq', column, value))
        return self
        
    def neq(self, column: str, value: Any) -> 'MockQuery':
        """Mock not equal filter."""
        self.filters.append(('neq', column, value))
        return self
        
    def gt(self, column: str, value: Any) -> 'MockQuery':
        """Mock greater than filter."""
        self.filters.append(('gt', column, value))
        return self
        
    def gte(self, column: str, value: Any) -> 'MockQuery':
        """Mock greater than or equal filter."""
        self.filters.append(('gte', column, value))
        return self
        
    def lt(self, column: str, value: Any) -> 'MockQuery':
        """Mock less than filter."""
        self.filters.append(('lt', column, value))
        return self
        
    def lte(self, column: str, value: Any) -> 'MockQuery':
        """Mock less than or equal filter."""
        self.filters.append(('lte', column, value))
        return self
        
    def like(self, column: str, pattern: str) -> 'MockQuery':
        """Mock LIKE filter."""
        self.filters.append(('like', column, pattern))
        return self

    def ilike(self, column: str, pattern: str) -> 'MockQuery':
        """Mock ILIKE filter."""
        self.filters.append(('ilike', column, pattern))
        return self

    def is_(self, column: str, value: Any) -> 'MockQuery':
        """Mock IS filter."""
        self.filters.append(('is', column, value))
        return self

    def in_(self, column: str, values: List) -> 'MockQuery':
        """Mock IN filter."""
        self.filters.append(('in', column, values))
        return self

    def contains(self, column: str, value: Any) -> 'MockQuery':
        """Mock contains filter for JSON columns."""
        self.filters.append(('contains', column, value))
        return self

    def order(self, column: str, desc: bool = False) -> 'MockQuery':
        """Mock order by."""
        self._order_by = (column, desc)
        return self

    def limit(self, count: int) -> 'MockQuery':
        """Mock limit."""
        self._limit = count
        return self

    def execute(self) -> MockQueryResult:
        """Mock query execution."""
        # Apply filters
        filtered_data = self.data.copy()
        
        for filter_type, column, value in self.filters:
            if filter_type == 'eq':
                filtered_data = [
                    item for item in filtered_data
                    if item.get(column) == value
                ]
            elif filter_type == 'neq':
                filtered_data = [
                    item for item in filtered_data
                    if item.get(column) != value
                ]
            elif filter_type == 'gt':
                filtered_data = [
                    item for item in filtered_data
                    if item.get(column, 0) > value
                ]
            elif filter_type == 'gte':
                filtered_data = [
                    item for item in filtered_data
                    if item.get(column, 0) >= value
                ]
            elif filter_type == 'lt':
                filtered_data = [
                    item for item in filtered_data
                    if item.get(column, 0) < value
                ]
            elif filter_type == 'lte':
                filtered_data = [
                    item for item in filtered_data
                    if item.get(column, 0) <= value
                ]
            elif filter_type == 'like':
                # Simple like implementation - just check if value is in column
                pattern = value.replace('%', '')
                filtered_data = [
                    item for item in filtered_data
                    if isinstance(item.get(column), str) and pattern in item.get(column, '')
                ]
            elif filter_type == 'ilike':
                # Case-insensitive like
                pattern = value.replace('%', '').lower()
                filtered_data = [
                    item for item in filtered_data
                    if isinstance(item.get(column), str) and pattern in item.get(column, '').lower()
                ]
            elif filter_type == 'is':
                filtered_data = [
                    item for item in filtered_data
                    if item.get(column) is value
                ]
            elif filter_type == 'in':
                filtered_data = [
                    item for item in filtered_data
                    if item.get(column) in value
                ]
            elif filter_type == 'contains':
                # For JSON columns that contain a value
                filtered_data = [
                    item for item in filtered_data
                    if column in item and value in item.get(column, [])
                ]
                
        # Apply order by
        if self._order_by:
            column, desc = self._order_by
            filtered_data = sorted(
                filtered_data,
                key=lambda x: (x.get(column) is None, x.get(column, '')),
                reverse=desc
            )
            
        # Apply limit
        if self._limit is not None and self._limit < len(filtered_data):
            filtered_data = filtered_data[:self._limit]
                
        return MockQueryResult(data=filtered_data)
        
    def insert(self, record_or_records: Union[Dict, List[Dict]]) -> 'MockQuery':
        """Mock insert method."""
        if isinstance(record_or_records, list):
            records = record_or_records
        else:
            records = [record_or_records]
            
        for record in records:
            # Add an ID if not present
            if 'id' not in record:
                import uuid
                record['id'] = str(uuid.uuid4())
                
            self.data.append(record.copy())
            
        return MockQuery(self.table, [records[-1]] if records else [])
        
    def update(self, record: Dict) -> 'MockQuery':
        """Mock update method."""
        updated_records = []
        
        for i, item in enumerate(self.data):
            match = True
            for filter_type, column, value in self.filters:
                if filter_type == 'eq' and item.get(column) != value:
                    match = False
                    break
                    
            if match:
                # Create updated record
                updated_record = {**item, **record}
                self.data[i] = updated_record
                updated_records.append(updated_record)
                
        return MockQuery(self.table, updated_records)
        
    def delete(self) -> 'MockQuery':
        """Mock delete method."""
        deleted_records = []
        remaining_records = []
        
        for item in self.data:
            should_delete = True
            for filter_type, column, value in self.filters:
                if filter_type == 'eq' and item.get(column) != value:
                    should_delete = False
                    break
                    
            if should_delete:
                deleted_records.append(item)
            else:
                remaining_records.append(item)
                
        self.data = remaining_records
        return MockQuery(self.table, deleted_records)

class MockStorage:
    """Mock Supabase Storage class."""
    
    def __init__(self):
        self.buckets = {}
        
    def from_(self, bucket: str) -> 'MockStorageBucket':
        """Get a bucket."""
        if bucket not in self.buckets:
            self.buckets[bucket] = {}
            
        return MockStorageBucket(bucket, self.buckets[bucket])

class MockStorageBucket:
    """Mock Supabase storage bucket."""
    
    def __init__(self, name: str, files: Dict[str, bytes]):
        self.name = name
        self.files = files
        
    def upload(self, path: str, file_contents: bytes) -> Dict:
        """Upload a file to the bucket."""
        self.files[path] = file_contents
        return {'Key': path}
        
    def download(self, path: str) -> bytes:
        """Download a file from the bucket."""
        if path not in self.files:
            raise ValueError(f"File not found: {path}")
            
        return self.files[path]
        
    def list(self, path: str = '') -> List[Dict]:
        """List files in the bucket."""
        result = []
        for file_path in self.files:
            if file_path.startswith(path):
                result.append({
                    'name': file_path,
                    'id': f"mock-id-{file_path}",
                    'updated_at': '2025-01-01T00:00:00Z'
                })
                
        return result
        
    def remove(self, paths: List[str]) -> Dict:
        """Remove files from the bucket."""
        for path in paths:
            if path in self.files:
                del self.files[path]
                
        return {'message': f"Deleted {len(paths)} files"}

class MockAuth:
    """Mock Supabase Auth class."""
    
    def __init__(self):
        self.users = {}
        self.current_user = None
        
    def sign_up(self, credentials: Dict) -> Dict:
        """Sign up a new user."""
        email = credentials.get('email')
        password = credentials.get('password')
        
        if not email or not password:
            return {'error': 'Email and password are required'}
            
        import uuid
        user_id = str(uuid.uuid4())
        
        user = {
            'id': user_id,
            'email': email,
            'user_metadata': credentials.get('options', {}).get('data', {})
        }
        
        self.users[user_id] = user
        self.current_user = user
        
        return {'user': user, 'session': {'access_token': f"mock-token-{user_id}"}}
        
    def sign_in(self, credentials: Dict) -> Dict:
        """Sign in an existing user."""
        email = credentials.get('email')
        password = credentials.get('password')
        
        if not email or not password:
            return {'error': 'Email and password are required'}
            
        # Find user by email
        user = next((u for u in self.users.values() if u['email'] == email), None)
        
        if not user:
            return {'error': 'Invalid credentials'}
            
        self.current_user = user
        
        return {'user': user, 'session': {'access_token': f"mock-token-{user['id']}"}}
        
    def user(self) -> Optional[Dict]:
        """Get the current user."""
        return self.current_user
        
    def session(self) -> Optional[Dict]:
        """Get the current session."""
        if not self.current_user:
            return None
            
        return {
            'access_token': f"mock-token-{self.current_user['id']}",
            'user': self.current_user
        }
        
    def update(self, attributes: Dict) -> Dict:
        """Update the current user."""
        if not self.current_user:
            return {'error': 'No active session'}
            
        user_id = self.current_user['id']
        self.users[user_id] = {**self.current_user, **attributes}
        self.current_user = self.users[user_id]
        
        return {'user': self.current_user}
        
    def sign_out(self) -> Dict:
        """Sign out the current user."""
        self.current_user = None
        return {'error': None}

class MockSupabase:
    """Mock Supabase client."""
    
    def __init__(self):
        self.tables = {}
        self.rpc_handlers = {}
        self.auth = MockAuth()
        self.storage = MockStorage()
        self.reset()
        
    def reset(self):
        """Reset mock state."""
        self.tables = {
            'molecules': [],
            'mixtures': [],
            'mixture_components': [],
            'migrations': [],
            'users': [],
            'experiments': [],
            'predictions': []
        }
        
        # Add default RPC handlers
        self.rpc_handlers = {
            'has_table': lambda params: {'data': [params.get('table_name') in self.tables]},
            'has_column': lambda params: {'data': [True]},  # Simplified - assume columns exist
            'exec_sql': lambda params: {'success': True}  # Default success for SQL execution
        }
        
    def table(self, name: str) -> MockQuery:
        """Get a query builder for a table."""
        if name not in self.tables:
            self.tables[name] = []
            
        return MockQuery(name, self.tables[name])
        
    def rpc(self, function: str, params: Optional[Dict] = None) -> MockQuery:
        """Mock RPC call."""
        params = params or {}
        
        if function in self.rpc_handlers:
            result = self.rpc_handlers[function](params)
            if isinstance(result, dict) and 'data' in result:
                return MockQuery('rpc', result['data'])
            return MockQuery('rpc', [result])
        else:
            return MockQuery('rpc', [{'success': True}])
            
    def sql(self, query: str, params: Optional[Dict] = None) -> MockQuery:
        """Mock SQL query."""
        # Simple query parser for testing
        query = query.lower()
        
        if query.startswith('select'):
            # Extract table name (very basic parser)
            parts = query.split('from')
            if len(parts) > 1:
                table_parts = parts[1].strip().split()
                if table_parts:
                    table_name = table_parts[0].strip().rstrip(';')
                    if table_name in self.tables:
                        return MockQuery(table_name, self.tables[table_name])
        
        # For modification queries, just return success
        return MockQuery('sql', [{'success': True}])
        
    def add_test_data(self, table: str, data: List[Dict]):
        """Add test data to a table."""
        if table not in self.tables:
            self.tables[table] = []
            
        self.tables[table].extend(data)
        
    def register_rpc_handler(self, function: str, handler: Callable):
        """Register a custom RPC handler."""
        self.rpc_handlers[function] = handler
        
    def from_(self, module: str):
        """Access different Supabase modules."""
        if module == 'storage':
            return self.storage
            
        raise ValueError(f"Unknown module: {module}")

@contextmanager
def patch_supabase_client(mock_client: Optional[MockSupabase] = None) -> Generator[MockSupabase, None, None]:
    """
    Context manager to patch the Supabase client.
    
    Args:
        mock_client: Optional mock Supabase client to use
        
    Yields:
        The mock client
    """
    if mock_client is None:
        mock_client = MockSupabase()
        
    # Targets to patch
    patch_targets = [
        'database.utils.connection.create_connection',
        'database.utils.connection.supabase.create_client',
        'supabase.create_client'
    ]
    
    patches = []
    for target in patch_targets:
        try:
            p = patch(target, return_value=mock_client)
            patches.append(p)
            p.start()
        except (ImportError, AttributeError):
            # Skip targets that don't exist
            pass
    
    try:
        yield mock_client
    finally:
        # Stop all patches
        for p in patches:
            p.stop()