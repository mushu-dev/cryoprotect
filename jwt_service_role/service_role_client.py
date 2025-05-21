"""
CryoProtect - Service Role Client

This module provides a client for interacting with the CryoProtect API using
service role authentication. It handles token acquisition, automatic refresh,
and request authentication.
"""

import time
import json
import logging
import threading
import requests
from typing import Dict, List, Optional, Union, Tuple, Any, Callable

from .token_manager import ServiceRoleTokenManager

# Configure logging
logger = logging.getLogger(__name__)

class ServiceRoleClient:
    """
    Client for making authenticated service role requests to the CryoProtect API.
    
    This client handles:
    - Acquiring and managing service role tokens
    - Automatically refreshing tokens when needed
    - Making authenticated API requests
    - Error handling and retries
    
    It includes thread-safety features to ensure tokens are properly managed
    in multi-threaded environments.
    """
    
    def __init__(self, 
                 base_url: str,
                 client_id: str,
                 client_scopes: List[str],
                 token_manager: Optional[ServiceRoleTokenManager] = None,
                 request_timeout: int = 30,
                 max_retries: int = 3,
                 retry_backoff: float = 1.5,
                 auto_refresh_threshold: int = 300):
        """
        Initialize the service role client.
        
        Args:
            base_url: Base URL for the API
            client_id: Client identifier for token requests
            client_scopes: List of scopes to request
            token_manager: Token manager instance (will create one if None)
            request_timeout: Timeout for API requests in seconds
            max_retries: Maximum number of retries for failed requests
            retry_backoff: Backoff factor for retries
            auto_refresh_threshold: Refresh token when this many seconds remain
        """
        self.base_url = base_url.rstrip('/')
        self.client_id = client_id
        self.client_scopes = client_scopes
        self.token_manager = token_manager or ServiceRoleTokenManager()
        self.request_timeout = request_timeout
        self.max_retries = max_retries
        self.retry_backoff = retry_backoff
        self.auto_refresh_threshold = auto_refresh_threshold
        
        # Token state
        self.current_token = None
        self.token_data = None
        self.token_expires_at = 0
        
        # Thread safety
        self.token_lock = threading.RLock()
        
        # Session for connection pooling
        self.session = requests.Session()
        
        # Initially acquire token
        self._acquire_token()
    
    def _acquire_token(self) -> bool:
        """
        Acquire a new service role token.
        
        Returns:
            True if token was acquired successfully, False otherwise
        """
        with self.token_lock:
            try:
                token, token_data = self.token_manager.generate_token(
                    client_id=self.client_id,
                    scopes=self.client_scopes
                )
                
                self.current_token = token
                self.token_data = token_data
                self.token_expires_at = token_data.get('exp', 0)
                
                logger.info(f"Acquired new service role token for client {self.client_id}")
                return True
            except Exception as e:
                logger.error(f"Failed to acquire service role token: {str(e)}")
                return False
    
    def _check_token_refresh(self) -> bool:
        """
        Check if token needs to be refreshed and refresh if needed.
        
        Returns:
            True if token is valid (refreshed or not), False otherwise
        """
        with self.token_lock:
            now = time.time()
            
            # If token is missing or expired, acquire a new one
            if not self.current_token or now >= self.token_expires_at:
                return self._acquire_token()
            
            # If token is approaching expiration, refresh it
            if self.token_expires_at - now < self.auto_refresh_threshold:
                logger.info("Token approaching expiration, refreshing")
                return self._acquire_token()
            
            # Token is valid and not approaching expiration
            return True
    
    def _get_auth_headers(self) -> Dict[str, str]:
        """
        Get authentication headers for API requests.
        
        Returns:
            Dictionary of headers including Authorization
        """
        with self.token_lock:
            return {
                'Authorization': f'Bearer {self.current_token}',
                'Content-Type': 'application/json'
            }
    
    def _make_request(self, 
                      method: str, 
                      endpoint: str, 
                      data: Optional[Dict] = None, 
                      params: Optional[Dict] = None, 
                      headers: Optional[Dict] = None,
                      files: Optional[Dict] = None,
                      stream: bool = False) -> requests.Response:
        """
        Make an authenticated API request with retries.
        
        Args:
            method: HTTP method (GET, POST, PUT, DELETE, etc.)
            endpoint: API endpoint path (will be joined with base_url)
            data: Request body data
            params: Query parameters
            headers: Additional headers
            files: Files to upload
            stream: Whether to stream the response
            
        Returns:
            Response object
            
        Raises:
            requests.RequestException: If the request fails after retries
        """
        # Ensure token is valid
        if not self._check_token_refresh():
            raise ValueError("Failed to acquire valid service role token")
        
        # Prepare request
        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        request_headers = self._get_auth_headers()
        if headers:
            request_headers.update(headers)
        
        # Convert data to JSON if it's a dict and not a file upload
        json_data = None
        if data and not files:
            json_data = data
            data = None
        
        # Make request with retries
        last_exception = None
        for attempt in range(self.max_retries + 1):
            try:
                response = self.session.request(
                    method=method,
                    url=url,
                    json=json_data,
                    data=data,
                    params=params,
                    headers=request_headers,
                    files=files,
                    timeout=self.request_timeout,
                    stream=stream
                )
                
                # Handle token expiration
                if response.status_code == 401:
                    logger.warning("Received 401, refreshing token")
                    self._acquire_token()
                    request_headers = self._get_auth_headers()
                    
                    # Retry immediately with new token if not the last attempt
                    if attempt < self.max_retries:
                        continue
                
                # Return response for other status codes
                return response
                
            except requests.RequestException as e:
                last_exception = e
                logger.warning(f"Request attempt {attempt + 1} failed: {str(e)}")
                
                # Break if this was the last attempt
                if attempt >= self.max_retries:
                    break
                
                # Backoff before retry
                backoff_time = self.retry_backoff ** attempt
                logger.info(f"Backing off for {backoff_time:.2f} seconds")
                time.sleep(backoff_time)
                
                # Refresh token if needed before retry
                self._check_token_refresh()
                request_headers = self._get_auth_headers()
        
        # All attempts failed
        if last_exception:
            raise last_exception
        
        # This shouldn't happen, but to satisfy type checking
        raise requests.RequestException("All request attempts failed")
    
    def get(self, 
            endpoint: str, 
            params: Optional[Dict] = None, 
            **kwargs) -> requests.Response:
        """
        Make a GET request to the API.
        
        Args:
            endpoint: API endpoint path
            params: Query parameters
            **kwargs: Additional arguments to pass to _make_request
            
        Returns:
            Response object
        """
        return self._make_request('GET', endpoint, params=params, **kwargs)
    
    def post(self, 
             endpoint: str, 
             data: Optional[Dict] = None, 
             **kwargs) -> requests.Response:
        """
        Make a POST request to the API.
        
        Args:
            endpoint: API endpoint path
            data: Request body data
            **kwargs: Additional arguments to pass to _make_request
            
        Returns:
            Response object
        """
        return self._make_request('POST', endpoint, data=data, **kwargs)
    
    def put(self, 
            endpoint: str, 
            data: Optional[Dict] = None, 
            **kwargs) -> requests.Response:
        """
        Make a PUT request to the API.
        
        Args:
            endpoint: API endpoint path
            data: Request body data
            **kwargs: Additional arguments to pass to _make_request
            
        Returns:
            Response object
        """
        return self._make_request('PUT', endpoint, data=data, **kwargs)
    
    def delete(self, 
               endpoint: str, 
               **kwargs) -> requests.Response:
        """
        Make a DELETE request to the API.
        
        Args:
            endpoint: API endpoint path
            **kwargs: Additional arguments to pass to _make_request
            
        Returns:
            Response object
        """
        return self._make_request('DELETE', endpoint, **kwargs)
    
    def patch(self, 
              endpoint: str, 
              data: Optional[Dict] = None, 
              **kwargs) -> requests.Response:
        """
        Make a PATCH request to the API.
        
        Args:
            endpoint: API endpoint path
            data: Request body data
            **kwargs: Additional arguments to pass to _make_request
            
        Returns:
            Response object
        """
        return self._make_request('PATCH', endpoint, data=data, **kwargs)
    
    def upload_file(self, 
                    endpoint: str, 
                    file_path: str, 
                    file_param_name: str = 'file',
                    additional_data: Optional[Dict] = None,
                    **kwargs) -> requests.Response:
        """
        Upload a file to the API.
        
        Args:
            endpoint: API endpoint path
            file_path: Path to the file to upload
            file_param_name: Name of the file parameter
            additional_data: Additional form data
            **kwargs: Additional arguments to pass to _make_request
            
        Returns:
            Response object
        """
        with open(file_path, 'rb') as f:
            files = {file_param_name: f}
            return self._make_request('POST', endpoint, data=additional_data, files=files, **kwargs)
    
    def download_file(self, 
                      endpoint: str, 
                      output_path: str,
                      params: Optional[Dict] = None,
                      **kwargs) -> bool:
        """
        Download a file from the API.
        
        Args:
            endpoint: API endpoint path
            output_path: Path to save the downloaded file
            params: Query parameters
            **kwargs: Additional arguments to pass to _make_request
            
        Returns:
            True if download was successful, False otherwise
        """
        try:
            response = self._make_request('GET', endpoint, params=params, stream=True, **kwargs)
            
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                return True
            else:
                logger.error(f"Download failed with status code {response.status_code}")
                return False
        except Exception as e:
            logger.error(f"Download failed: {str(e)}")
            return False
    
    def close(self):
        """
        Close the client session.
        """
        self.session.close()
    
    def __enter__(self):
        """
        Enter context manager.
        """
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit context manager.
        """
        self.close()