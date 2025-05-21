#!/usr/bin/env python3
"""
Service Role Client for CryoProtect.

This module provides a client for making authenticated API requests using
service role authentication with JWT tokens.
"""

import os
import time
import json
import logging
import requests
from typing import Dict, Any, Optional, Union, List
from urllib.parse import urljoin

from .token_manager import create_service_token, validate_service_token

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ServiceRoleClient:
    """
    Client for making authenticated API requests using service role authentication.
    
    Features:
    - Automatic token acquisition and renewal
    - Request retries with exponential backoff
    - Connection pooling for better performance
    - Request/response logging
    """
    
    def __init__(
        self,
        base_url: str,
        service_name: str,
        scopes: Optional[List[str]] = None,
        token_expiration: int = 3600,
        max_retries: int = 3,
        timeout: int = 30,
        session: Optional[requests.Session] = None,
        verify_ssl: bool = True
    ):
        """
        Initialize the service role client.
        
        Args:
            base_url: Base URL for API requests
            service_name: Name of the service making requests
            scopes: Permission scopes for the service
            token_expiration: Token expiration time in seconds
            max_retries: Maximum number of request retries
            timeout: Request timeout in seconds
            session: Optional requests.Session for connection pooling
            verify_ssl: Whether to verify SSL certificates
        """
        self.base_url = base_url.rstrip('/') + '/'
        self.service_name = service_name
        self.scopes = scopes or ['*']
        self.token_expiration = token_expiration
        self.max_retries = max_retries
        self.timeout = timeout
        self.verify_ssl = verify_ssl
        
        # Use provided session or create a new one
        self.session = session or requests.Session()
        
        # Token cache
        self.current_token = None
        self.token_expiry = 0
        
        # Get initial token
        self._refresh_token()
    
    def _refresh_token(self):
        """Get a fresh service token."""
        try:
            self.current_token = create_service_token(
                service_name=self.service_name,
                scopes=self.scopes,
                expiration=self.token_expiration
            )
            
            # Set token expiry to slightly before actual expiration
            self.token_expiry = time.time() + (self.token_expiration * 0.9)
            
            logger.debug(f"Acquired new service token for {self.service_name}")
        except Exception as e:
            logger.error(f"Error refreshing service token: {str(e)}")
            raise
    
    def _ensure_valid_token(self):
        """Ensure the current token is valid, refreshing if needed."""
        if not self.current_token or time.time() >= self.token_expiry:
            self._refresh_token()
    
    def _get_auth_headers(self) -> Dict[str, str]:
        """Get authentication headers with the current token."""
        self._ensure_valid_token()
        return {
            'Authorization': f'Bearer {self.current_token}'
        }
    
    def request(
        self,
        method: str,
        path: str,
        params: Optional[Dict[str, Any]] = None,
        data: Optional[Any] = None,
        json_data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None,
        retry_count: int = 0
    ) -> requests.Response:
        """
        Make an authenticated API request.
        
        Args:
            method: HTTP method (GET, POST, PUT, DELETE, etc.)
            path: API endpoint path (relative to base_url)
            params: Query parameters
            data: Request body data
            json_data: JSON request body
            headers: Additional request headers
            timeout: Request timeout in seconds (overrides default)
            retry_count: Current retry count (used internally)
            
        Returns:
            requests.Response object
            
        Raises:
            requests.RequestException: If the request fails
        """
        # Build URL
        url = urljoin(self.base_url, path.lstrip('/'))
        
        # Prepare headers
        request_headers = headers.copy() if headers else {}
        request_headers.update(self._get_auth_headers())
        
        # Set content type for JSON data
        if json_data and 'Content-Type' not in request_headers:
            request_headers['Content-Type'] = 'application/json'
        
        # Set default timeout
        request_timeout = timeout if timeout is not None else self.timeout
        
        try:
            # Make the request
            start_time = time.time()
            
            response = self.session.request(
                method=method,
                url=url,
                params=params,
                data=data,
                json=json_data,
                headers=request_headers,
                timeout=request_timeout,
                verify=self.verify_ssl
            )
            
            elapsed = time.time() - start_time
            
            # Log request details
            logger.debug(
                f"{method} {url} - Status: {response.status_code}, "
                f"Time: {elapsed:.3f}s, Size: {len(response.content)} bytes"
            )
            
            # Check if we need to retry (server errors)
            if response.status_code >= 500 and retry_count < self.max_retries:
                return self._retry_request(
                    method, path, params, data, json_data,
                    headers, timeout, retry_count + 1
                )
            
            # Check for authentication errors (refresh token and retry)
            if response.status_code == 401 and retry_count < self.max_retries:
                self._refresh_token()
                return self._retry_request(
                    method, path, params, data, json_data,
                    headers, timeout, retry_count + 1
                )
            
            # Raise for status
            response.raise_for_status()
            
            return response
        
        except requests.RequestException as e:
            # Retry on connection errors
            if retry_count < self.max_retries and isinstance(
                e, (requests.ConnectionError, requests.Timeout)
            ):
                return self._retry_request(
                    method, path, params, data, json_data,
                    headers, timeout, retry_count + 1
                )
            
            logger.error(f"Request failed: {str(e)}")
            raise
    
    def _retry_request(
        self,
        method: str,
        path: str,
        params: Optional[Dict[str, Any]],
        data: Optional[Any],
        json_data: Optional[Dict[str, Any]],
        headers: Optional[Dict[str, str]],
        timeout: Optional[int],
        retry_count: int
    ) -> requests.Response:
        """Retry a failed request with exponential backoff."""
        # Calculate delay with exponential backoff
        delay = 0.5 * (2 ** (retry_count - 1))
        
        logger.warning(
            f"Retrying {method} {path} (attempt {retry_count}/{self.max_retries}) "
            f"after {delay:.2f}s delay"
        )
        
        time.sleep(delay)
        
        return self.request(
            method, path, params, data, json_data,
            headers, timeout, retry_count
        )
    
    def get(
        self,
        path: str,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None
    ) -> requests.Response:
        """
        Make a GET request.
        
        Args:
            path: API endpoint path
            params: Query parameters
            headers: Additional request headers
            timeout: Request timeout in seconds
            
        Returns:
            requests.Response object
        """
        return self.request(
            method='GET',
            path=path,
            params=params,
            headers=headers,
            timeout=timeout
        )
    
    def post(
        self,
        path: str,
        data: Optional[Any] = None,
        json_data: Optional[Dict[str, Any]] = None,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None
    ) -> requests.Response:
        """
        Make a POST request.
        
        Args:
            path: API endpoint path
            data: Request body data
            json_data: JSON request body
            params: Query parameters
            headers: Additional request headers
            timeout: Request timeout in seconds
            
        Returns:
            requests.Response object
        """
        return self.request(
            method='POST',
            path=path,
            params=params,
            data=data,
            json_data=json_data,
            headers=headers,
            timeout=timeout
        )
    
    def put(
        self,
        path: str,
        data: Optional[Any] = None,
        json_data: Optional[Dict[str, Any]] = None,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None
    ) -> requests.Response:
        """
        Make a PUT request.
        
        Args:
            path: API endpoint path
            data: Request body data
            json_data: JSON request body
            params: Query parameters
            headers: Additional request headers
            timeout: Request timeout in seconds
            
        Returns:
            requests.Response object
        """
        return self.request(
            method='PUT',
            path=path,
            params=params,
            data=data,
            json_data=json_data,
            headers=headers,
            timeout=timeout
        )
    
    def patch(
        self,
        path: str,
        data: Optional[Any] = None,
        json_data: Optional[Dict[str, Any]] = None,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None
    ) -> requests.Response:
        """
        Make a PATCH request.
        
        Args:
            path: API endpoint path
            data: Request body data
            json_data: JSON request body
            params: Query parameters
            headers: Additional request headers
            timeout: Request timeout in seconds
            
        Returns:
            requests.Response object
        """
        return self.request(
            method='PATCH',
            path=path,
            params=params,
            data=data,
            json_data=json_data,
            headers=headers,
            timeout=timeout
        )
    
    def delete(
        self,
        path: str,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: Optional[int] = None
    ) -> requests.Response:
        """
        Make a DELETE request.
        
        Args:
            path: API endpoint path
            params: Query parameters
            headers: Additional request headers
            timeout: Request timeout in seconds
            
        Returns:
            requests.Response object
        """
        return self.request(
            method='DELETE',
            path=path,
            params=params,
            headers=headers,
            timeout=timeout
        )

# Create default client factory
def create_client(
    service_name: str,
    base_url: Optional[str] = None,
    scopes: Optional[List[str]] = None
) -> ServiceRoleClient:
    """
    Create a service role client.
    
    Args:
        service_name: Name of the service making requests
        base_url: Base URL for API requests (defaults to environment variable)
        scopes: Permission scopes for the service
        
    Returns:
        ServiceRoleClient instance
    """
    # Get base URL from environment if not provided
    url = base_url or os.environ.get('API_BASE_URL')
    
    if not url:
        raise ValueError("Base URL must be provided or set in API_BASE_URL environment variable")
    
    return ServiceRoleClient(
        base_url=url,
        service_name=service_name,
        scopes=scopes
    )

# Main function for CLI usage
def main():
    """Main function for command-line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Service Role Client')
    parser.add_argument('--service', required=True, help='Service name')
    parser.add_argument('--url', required=True, help='API base URL')
    parser.add_argument('--method', default='GET', help='HTTP method')
    parser.add_argument('--path', required=True, help='API endpoint path')
    parser.add_argument('--params', help='Query parameters (JSON)')
    parser.add_argument('--data', help='Request body data')
    parser.add_argument('--json', help='JSON request body (JSON string)')
    parser.add_argument('--headers', help='Additional headers (JSON)')
    parser.add_argument('--timeout', type=int, default=30, help='Request timeout')
    parser.add_argument('--retries', type=int, default=3, help='Max retries')
    
    args = parser.parse_args()
    
    # Create client
    client = ServiceRoleClient(
        base_url=args.url,
        service_name=args.service,
        max_retries=args.retries,
        timeout=args.timeout
    )
    
    # Parse JSON arguments
    params = json.loads(args.params) if args.params else None
    json_data = json.loads(args.json) if args.json else None
    headers = json.loads(args.headers) if args.headers else None
    
    # Make request
    response = client.request(
        method=args.method.upper(),
        path=args.path,
        params=params,
        data=args.data,
        json_data=json_data,
        headers=headers,
        timeout=args.timeout
    )
    
    # Print response
    print(f"Status Code: {response.status_code}")
    print(f"Headers: {json.dumps(dict(response.headers), indent=2)}")
    
    # Try to parse JSON response
    try:
        response_data = response.json()
        print(f"Response: {json.dumps(response_data, indent=2)}")
    except ValueError:
        print(f"Response: {response.text}")

if __name__ == "__main__":
    main()