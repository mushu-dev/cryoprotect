"""
CryoProtect Analyzer - Mock Supabase Query

This module provides mock implementations of Supabase query builders.
"""

import uuid
from datetime import datetime
from .data import _mock_data, _mock_views, _update_molecule_with_properties_view

# Mock Response class
class MockResponse:
    """Mock response object that mimics Supabase response."""
    
    def __init__(self, data=None, error=None, count=None):
        self.data = data if data is not None else []
        self.error = error
        self.count = count
        
        # For auth responses
        self.user = None
        self.session = None
        
        if data and isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
            if 'email' in data[0]:
                self.user = data[0]
    
    def __repr__(self):
        return f"MockResponse(data={self.data}, error={self.error}, count={self.count})"


# Mock Query Builder
class MockQueryBuilder:
    """Mock query builder that mimics Supabase query builder."""
    
    def __init__(self, table_name):
        self.table_name = table_name
        self.select_columns = '*'
        self.filters = []
        self.range_start = None
        self.range_end = None
        self.limit_val = None
        self.order_columns = []
        self.count_requested = False
    
    def select(self, *columns, count=None):
        if columns and columns[0] != '*':
            self.select_columns = columns
        if count:
            self.count_requested = True
        return self
    
    def eq(self, column, value):
        self.filters.append(('eq', column, value))
        return self
    
    def neq(self, column, value):
        self.filters.append(('neq', column, value))
        return self
    
    def gt(self, column, value):
        self.filters.append(('gt', column, value))
        return self
    
    def gte(self, column, value):
        self.filters.append(('gte', column, value))
        return self
    
    def lt(self, column, value):
        self.filters.append(('lt', column, value))
        return self
    
    def lte(self, column, value):
        self.filters.append(('lte', column, value))
        return self
    
    def like(self, column, value):
        self.filters.append(('like', column, value))
        return self
    
    def ilike(self, column, value):
        self.filters.append(('ilike', column, value))
        return self
    
    def is_(self, column, value):
        self.filters.append(('is', column, value))
        return self
    
    def in_(self, column, values):
        self.filters.append(('in', column, values))
        return self
    
    def range(self, start, end):
        self.range_start = start
        self.range_end = end
        return self
    
    def limit(self, count):
        self.limit_val = count
        return self
    
    def order(self, column, ascending=True):
        self.order_columns.append((column, ascending))
        return self
    
    def execute(self):
        """Execute the query and return a mock response."""
        try:
            # Get data from the table
            if self.table_name in _mock_data:
                data = _mock_data[self.table_name]
            elif self.table_name in _mock_views:
                data = _mock_views[self.table_name]
            else:
                return MockResponse([], None)
            
            # Apply filters
            filtered_data = self._apply_filters(data)
            
            # Apply ordering
            if self.order_columns:
                for column, ascending in reversed(self.order_columns):
                    filtered_data = sorted(
                        filtered_data,
                        key=lambda x: x.get(column, ''),
                        reverse=not ascending
                    )
            
            # Apply range/pagination
            if self.range_start is not None and self.range_end is not None:
                filtered_data = filtered_data[self.range_start:self.range_end + 1]
            elif self.limit_val:
                filtered_data = filtered_data[:self.limit_val]
            
            # Handle count
            count = len(filtered_data) if self.count_requested else None
            
            # Handle select columns
            if self.select_columns != '*':
                filtered_data = [
                    {col: item.get(col) for col in self.select_columns}
                    for item in filtered_data
                ]
            
            return MockResponse(filtered_data, None, count)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})
    
    def _apply_filters(self, data):
        """Apply filters to the data."""
        result = data
        
        for filter_type, column, value in self.filters:
            if filter_type == 'eq':
                result = [item for item in result if item.get(column) == value]
            elif filter_type == 'neq':
                result = [item for item in result if item.get(column) != value]
            elif filter_type == 'gt':
                result = [item for item in result if item.get(column, 0) > value]
            elif filter_type == 'gte':
                result = [item for item in result if item.get(column, 0) >= value]
            elif filter_type == 'lt':
                result = [item for item in result if item.get(column, 0) < value]
            elif filter_type == 'lte':
                result = [item for item in result if item.get(column, 0) <= value]
            elif filter_type == 'like':
                if isinstance(value, str):
                    value = value.replace('%', '')
                    result = [item for item in result if value in str(item.get(column, ''))]
            elif filter_type == 'ilike':
                if isinstance(value, str):
                    value = value.replace('%', '').lower()
                    result = [item for item in result if value in str(item.get(column, '')).lower()]
            elif filter_type == 'is':
                result = [item for item in result if item.get(column) is value]
            elif filter_type == 'in':
                result = [item for item in result if item.get(column) in value]
        
        return result
    
    def insert(self, data):
        """Prepare to insert data."""
        return MockInsertBuilder(self.table_name, data)
    
    def update(self, data):
        """Prepare to update data."""
        return MockUpdateBuilder(self.table_name, data, self.filters)
    
    def upsert(self, data):
        """Prepare to upsert data."""
        return MockUpsertBuilder(self.table_name, data)
    
    def delete(self):
        """Prepare to delete data."""
        return MockDeleteBuilder(self.table_name, self.filters)


class MockInsertBuilder:
    """Mock insert builder."""
    
    def __init__(self, table_name, data):
        self.table_name = table_name
        self.data = data
    
    def execute(self):
        """Execute the insert operation."""
        try:
            if self.table_name not in _mock_data:
                _mock_data[self.table_name] = []
            
            # Handle single item or list
            items_to_insert = self.data if isinstance(self.data, list) else [self.data]
            inserted_items = []
            
            for item in items_to_insert:
                # Ensure item has an ID
                if 'id' not in item:
                    item['id'] = str(uuid.uuid4())
                
                # Add timestamps if not present
                now = datetime.now().isoformat()
                if 'created_at' not in item:
                    item['created_at'] = now
                if 'updated_at' not in item:
                    item['updated_at'] = now
                
                # Add to the table
                _mock_data[self.table_name].append(item)
                inserted_items.append(item)
            
            # Update views if needed
            if self.table_name in ['molecules', 'molecular_properties']:
                _update_molecule_with_properties_view()
            
            return MockResponse(inserted_items)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})


class MockUpdateBuilder:
    """Mock update builder."""
    
    def __init__(self, table_name, data, filters):
        self.table_name = table_name
        self.data = data
        self.filters = filters
    
    def eq(self, column, value):
        self.filters.append(('eq', column, value))
        return self
    
    def execute(self):
        """Execute the update operation."""
        try:
            if self.table_name not in _mock_data:
                return MockResponse([], {"message": f"Table {self.table_name} not found"})
            
            # Find items to update
            query_builder = MockQueryBuilder(self.table_name)
            query_builder.filters = self.filters
            response = query_builder.execute()
            
            if response.error:
                return response
            
            items_to_update = response.data
            updated_items = []
            
            for item in items_to_update:
                # Find the item in the original data
                original_item = next(
                    (i for i in _mock_data[self.table_name] if i['id'] == item['id']),
                    None
                )
                
                if original_item:
                    # Update the item
                    for key, value in self.data.items():
                        original_item[key] = value
                    
                    # Update timestamp
                    original_item['updated_at'] = datetime.now().isoformat()
                    
                    updated_items.append(original_item)
            
            # Update views if needed
            if self.table_name in ['molecules', 'molecular_properties']:
                _update_molecule_with_properties_view()
            
            return MockResponse(updated_items)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})


class MockUpsertBuilder:
    """Mock upsert builder."""
    
    def __init__(self, table_name, data):
        self.table_name = table_name
        self.data = data
    
    def execute(self):
        """Execute the upsert operation."""
        try:
            if self.table_name not in _mock_data:
                _mock_data[self.table_name] = []
            
            # Handle single item or list
            items_to_upsert = self.data if isinstance(self.data, list) else [self.data]
            upserted_items = []
            
            for item in items_to_upsert:
                # Check if item exists
                existing_item = None
                if 'id' in item:
                    existing_item = next(
                        (i for i in _mock_data[self.table_name] if i['id'] == item['id']),
                        None
                    )
                
                if existing_item:
                    # Update existing item
                    for key, value in item.items():
                        existing_item[key] = value
                    
                    # Update timestamp
                    existing_item['updated_at'] = datetime.now().isoformat()
                    
                    upserted_items.append(existing_item)
                else:
                    # Insert new item
                    if 'id' not in item:
                        item['id'] = str(uuid.uuid4())
                    
                    # Add timestamps if not present
                    now = datetime.now().isoformat()
                    if 'created_at' not in item:
                        item['created_at'] = now
                    if 'updated_at' not in item:
                        item['updated_at'] = now
                    
                    # Add to the table
                    _mock_data[self.table_name].append(item)
                    upserted_items.append(item)
            
            # Update views if needed
            if self.table_name in ['molecules', 'molecular_properties']:
                _update_molecule_with_properties_view()
            
            return MockResponse(upserted_items)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})


class MockDeleteBuilder:
    """Mock delete builder."""
    
    def __init__(self, table_name, filters):
        self.table_name = table_name
        self.filters = filters
    
    def eq(self, column, value):
        self.filters.append(('eq', column, value))
        return self
    
    def execute(self):
        """Execute the delete operation."""
        try:
            if self.table_name not in _mock_data:
                return MockResponse([], {"message": f"Table {self.table_name} not found"})
            
            # Find items to delete
            query_builder = MockQueryBuilder(self.table_name)
            query_builder.filters = self.filters
            response = query_builder.execute()
            
            if response.error:
                return response
            
            items_to_delete = response.data
            deleted_ids = [item['id'] for item in items_to_delete]
            
            # Remove items from the table
            _mock_data[self.table_name] = [
                item for item in _mock_data[self.table_name] 
                if item['id'] not in deleted_ids
            ]
            
            # Update views if needed
            if self.table_name in ['molecules', 'molecular_properties']:
                _update_molecule_with_properties_view()
            
            return MockResponse(items_to_delete)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})