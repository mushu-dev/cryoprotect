"""
Convex database adapter for CryoProtect.

This module provides a compatibility layer between the existing Supabase-based
API and the new Convex database backend.
"""

import os
import json
import requests
import logging
from urllib.parse import urljoin

logger = logging.getLogger(__name__)

class ConvexAdapter:
    """
    Adapter class that provides Supabase-like interface for Convex.
    This allows for gradual migration of the codebase.
    """
    
    def __init__(self, url=None, key=None):
        """
        Initialize the Convex adapter.
        
        Args:
            url (str): The Convex URL. If not provided, it's read from CONVEX_URL env var.
            key (str): The Convex deployment key. If not provided, it's read from CONVEX_DEPLOYMENT_KEY env var.
        """
        self.url = url or os.environ.get('CONVEX_URL', 'https://dynamic-mink-63.convex.cloud')
        self.key = key or os.environ.get('CONVEX_DEPLOYMENT_KEY', '')
        
        # Ensure the URL ends with a slash for proper joining
        if not self.url.endswith('/'):
            self.url = self.url + '/'
            
        self.headers = {
            'Content-Type': 'application/json',
            'Authorization': f'Bearer {self.key}' if self.key else None
        }
        
        logger.info("Initialized Convex adapter with URL: %s", self.url)
    
    def table(self, table_name):
        """
        Get a reference to a Convex table (collection).
        
        Args:
            table_name (str): The name of the table.
            
        Returns:
            TableAdapter: An adapter for the specified table.
        """
        return TableAdapter(self, table_name)
    
    def auth(self):
        """
        Get a reference to the Convex auth methods.
        
        Returns:
            AuthAdapter: An adapter for authentication methods.
        """
        return AuthAdapter(self)
    
    def execute_query(self, action, path, params=None):
        """
        Execute a query against Convex.
        
        Args:
            action (str): The HTTP method (GET, POST, PUT, DELETE).
            path (str): The API path.
            params (dict): The parameters to send with the request.
            
        Returns:
            dict: The response from Convex.
        """
        url = urljoin(self.url, path)
        
        try:
            if action == 'GET':
                response = requests.get(url, params=params, headers=self.headers)
            elif action == 'POST':
                response = requests.post(url, json=params, headers=self.headers)
            elif action == 'PUT':
                response = requests.put(url, json=params, headers=self.headers)
            elif action == 'DELETE':
                response = requests.delete(url, json=params, headers=self.headers)
            else:
                raise ValueError(f"Unsupported action: {action}")
            
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error("Error executing Convex query: %s", str(e))
            raise

class TableAdapter:
    """Adapter for Convex tables to provide Supabase-like interface."""
    
    def __init__(self, adapter, table_name):
        """
        Initialize the table adapter.
        
        Args:
            adapter (ConvexAdapter): The parent Convex adapter.
            table_name (str): The name of the table.
        """
        self.adapter = adapter
        self.table_name = table_name
        self.filters = {}
        self.order_clauses = []
        self.limit_value = None
        self.offset_value = None
        
    def select(self, columns="*"):
        """
        Select columns from the table.
        
        Args:
            columns (str): Comma-separated column names or "*" for all columns.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.columns = columns
        return self
    
    def eq(self, column, value):
        """
        Add an equality filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[column] = value
        return self
    
    def neq(self, column, value):
        """
        Add a not-equal filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_neq"] = value
        return self
    
    def gt(self, column, value):
        """
        Add a greater-than filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_gt"] = value
        return self
    
    def gte(self, column, value):
        """
        Add a greater-than-or-equal filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_gte"] = value
        return self
    
    def lt(self, column, value):
        """
        Add a less-than filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_lt"] = value
        return self
    
    def lte(self, column, value):
        """
        Add a less-than-or-equal filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_lte"] = value
        return self
    
    def order(self, column, order="asc"):
        """
        Add an ordering clause.
        
        Args:
            column (str): The column name.
            order (str): Either "asc" or "desc".
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.order_clauses.append((column, order))
        return self
    
    def limit(self, limit):
        """
        Set a limit on the number of results.
        
        Args:
            limit (int): The maximum number of results to return.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.limit_value = limit
        return self
    
    def offset(self, offset):
        """
        Set an offset for the results.
        
        Args:
            offset (int): The offset for results.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.offset_value = offset
        return self
    
    def execute(self):
        """
        Execute the query.
        
        Returns:
            Response: A Supabase-like response object.
        """
        # Construct query parameters
        params = {
            "table": self.table_name,
            "filters": self.filters,
        }
        
        if hasattr(self, 'columns'):
            params["columns"] = self.columns
            
        if self.order_clauses:
            params["order"] = self.order_clauses
            
        if self.limit_value is not None:
            params["limit"] = self.limit_value
            
        if self.offset_value is not None:
            params["offset"] = self.offset_value
        
        # Execute the query and convert to Supabase-like response
        try:
            result = self.adapter.execute_query('POST', 'api/query', params)
            return Response(result.get('data', []), result.get('error'))
        except Exception as e:
            logger.error("Error executing Convex query: %s", str(e))
            return Response([], str(e))
    
    def insert(self, data):
        """
        Insert data into the table.
        
        Args:
            data (dict or list): The data to insert.
            
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "table": self.table_name,
            "data": data
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/insert', params)
            return Response(result.get('data', []), result.get('error'))
        except Exception as e:
            logger.error("Error inserting into Convex: %s", str(e))
            return Response([], str(e))
    
    def update(self, data):
        """
        Update data in the table.
        
        Args:
            data (dict): The data to update.
            
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "table": self.table_name,
            "data": data,
            "filters": self.filters
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/update', params)
            return Response(result.get('data', []), result.get('error'))
        except Exception as e:
            logger.error("Error updating in Convex: %s", str(e))
            return Response([], str(e))
    
    def delete(self):
        """
        Delete data from the table.
        
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "table": self.table_name,
            "filters": self.filters
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/delete', params)
            return Response(result.get('data', []), result.get('error'))
        except Exception as e:
            logger.error("Error deleting from Convex: %s", str(e))
            return Response([], str(e))

class AuthAdapter:
    """Adapter for Convex auth to provide Supabase-like interface."""
    
    def __init__(self, adapter):
        """
        Initialize the auth adapter.
        
        Args:
            adapter (ConvexAdapter): The parent Convex adapter.
        """
        self.adapter = adapter
    
    def sign_in_with_password(self, credentials):
        """
        Sign in with email and password.
        
        Args:
            credentials (dict): Dict with email and password.
            
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "email": credentials.get('email'),
            "password": credentials.get('password')
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/auth/signin', params)
            return Response(result.get('data', {}), result.get('error'))
        except Exception as e:
            logger.error("Error signing in with Convex: %s", str(e))
            return Response({}, str(e))
    
    def sign_up(self, credentials):
        """
        Sign up with email and password.
        
        Args:
            credentials (dict): Dict with email and password.
            
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "email": credentials.get('email'),
            "password": credentials.get('password')
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/auth/signup', params)
            return Response(result.get('data', {}), result.get('error'))
        except Exception as e:
            logger.error("Error signing up with Convex: %s", str(e))
            return Response({}, str(e))
    
    def sign_out(self):
        """
        Sign out the current user.
        
        Returns:
            Response: A Supabase-like response object.
        """
        try:
            result = self.adapter.execute_query('POST', 'api/auth/signout', {})
            return Response(result.get('data', {}), result.get('error'))
        except Exception as e:
            logger.error("Error signing out with Convex: %s", str(e))
            return Response({}, str(e))

class Response:
    """A Supabase-like response object."""
    
    def __init__(self, data, error=None):
        """
        Initialize the response.
        
        Args:
            data: The response data.
            error: Any error that occurred.
        """
        self.data = data
        self.error = error

def create_client(url=None, key=None, use_convex=None):
    """
    Create a Convex client based on environment configuration.
    
    Args:
        url (str): The Convex URL.
        key (str): The Convex deployment key.
        use_convex (bool): Force use of Convex if True, Supabase if False.
        
    Returns:
        ConvexAdapter or supabase.Client: The appropriate client.
    """
    # Determine whether to use Convex or fall back to Supabase
    use_convex_env = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
    should_use_convex = use_convex if use_convex is not None else use_convex_env
    
    if should_use_convex:
        logger.info("Using Convex adapter for database operations")
        return ConvexAdapter(url, key)
    else:
        # Fall back to Supabase
        logger.info("Using Supabase for database operations")
        from supabase import create_client as create_supabase_client
        
        supabase_url = os.environ.get('SUPABASE_URL', '')
        supabase_key = os.environ.get('SUPABASE_KEY', '')
        
        return create_supabase_client(supabase_url, supabase_key)