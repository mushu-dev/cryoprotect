# Supabase Connection Implementation Guide

This document provides a comprehensive guide for implementing a robust solution to the database connection issues in the CryoProtect v2 project. It's designed for Roo code agents to understand the full context and implement a solution that addresses all aspects of the problem.

## 1. Problem Analysis

### 1.1 Current Issues

Based on logs and code analysis, we have identified these specific issues:

1. **DNS Resolution Failure**: The application cannot resolve the hostname `db.tsdlmynydfuypiugmkev.supabase.co`
   - Error: `could not translate host name "db.tsdlmynydfuypiugmkev.supabase.co" to address: No such host is known`

2. **Configuration Inconsistency**: Different environment variable patterns (`DB_*` vs `SUPABASE_DB_*`)
   - Some files use `DB_HOST`, others use `SUPABASE_DB_HOST`
   - The `get_db_config()` function only looks for `DB_*` variables

3. **Authentication Failure**: When falling back to local database, authentication fails
   - Error: `FATAL: password authentication failed for user "postgres"`

4. **Connection Pool Issues**: The current connection pool wrapper doesn't handle DNS failures gracefully

### 1.2 Project Environment

- **Project ID**: `tsdlmynydfuypiugmkev`
- **Database Host**: `db.tsdlmynydfuypiugmkev.supabase.co`
- **Project Status**: Active and healthy according to Supabase MCP
- **Required Files**: connection_pool_wrapper.py, config.py, .env
- **Verification Script**: verify_imported_data.py

## 2. Solution Architecture

The solution consists of four key components:

1. **IP Resolution**: Multiple methods to resolve the database hostname to an IP address
2. **Configuration Unification**: Ensure consistent environment variables
3. **Enhanced Connection Pool**: A robust connection wrapper that handles DNS failures
4. **Testing & Verification**: Comprehensive testing to validate fixes

### 2.1 Component Relationships

```
┌─────────────────┐        ┌──────────────────┐        ┌────────────────┐
│  IP Resolution  │───────▶│  Config Updates  │───────▶│  Connection    │
└─────────────────┘        └──────────────────┘        │  Pool Wrapper  │
        │                          │                   └────────────────┘
        │                          │                          │
        │                          │                          │
        ▼                          ▼                          ▼
┌──────────────────────────────────────────────────────────────────────┐
│                     Verification & Testing                            │
└──────────────────────────────────────────────────────────────────────┘
```

## 3. Implementation Guide

### 3.1 IP Resolution Implementation

#### 3.1.1 Methods to Implement

1. **Standard DNS Resolution**: Use socket.gethostbyname()
2. **Alternative DNS Servers**: Try multiple DNS servers (8.8.8.8, 1.1.1.1)
3. **MCP Resolution**: Use Supabase MCP to get project info
4. **Heuristic Method**: Use patterns in IP assignment as fallback

#### 3.1.2 Code Example: Multi-Method Resolution

```python
def resolve_ip_address(hostname: str) -> Optional[str]:
    """Attempt to resolve IP address using multiple methods."""
    # Method 1: Standard resolution
    try:
        ip_address = socket.gethostbyname(hostname)
        logger.info(f"Standard DNS resolution succeeded: {ip_address}")
        return ip_address
    except socket.gaierror:
        logger.warning("Standard DNS resolution failed")
    
    # Method 2: Alternative DNS servers
    for dns_server in ["8.8.8.8", "1.1.1.1"]:
        try:
            # Windows vs Unix implementation differences here
            # On Windows: Use nslookup
            if platform.system() == "Windows":
                result = subprocess.run(
                    ["nslookup", hostname, dns_server],
                    capture_output=True, text=True
                )
                # Parse output...
            else:
                # On Unix: Use dig
                result = subprocess.run(
                    ["dig", "+short", hostname, f"@{dns_server}"],
                    capture_output=True, text=True
                )
                if result.stdout.strip():
                    ip_address = result.stdout.strip().split("\n")[0]
                    return ip_address
        except Exception:
            continue
    
    # Method 3: MCP Resolution
    # Implementation details...
    
    # Method 4: Heuristic fallback
    # Implementation details...
```

### 3.2 Configuration Unification

#### 3.2.1 Environment Variables to Set

- `SUPABASE_URL`
- `SUPABASE_PROJECT_ID`
- `SUPABASE_DB_HOST`
- `SUPABASE_DB_PORT`
- `SUPABASE_DB_NAME`
- `SUPABASE_DB_USER`
- `SUPABASE_DB_PASSWORD`
- `DB_HOST`
- `DB_PORT`
- `DB_NAME`
- `DB_USER`
- `DB_PASSWORD`
- `SUPABASE_DB_IP` (new variable for direct IP connection)

#### 3.2.2 Code Example: Loading Environment Variables

```python
def load_environment_variables() -> Dict[str, str]:
    """Load and normalize environment variables."""
    load_dotenv()
    
    env_vars = {
        # Supabase connection variables
        "SUPABASE_URL": os.getenv("SUPABASE_URL", ""),
        "SUPABASE_KEY": os.getenv("SUPABASE_KEY", ""),
        "SUPABASE_PROJECT_ID": os.getenv("SUPABASE_PROJECT_ID", ""),
        
        # Database connection variables (SUPABASE_DB_ prefix)
        "SUPABASE_DB_HOST": os.getenv("SUPABASE_DB_HOST", ""),
        "SUPABASE_DB_PORT": os.getenv("SUPABASE_DB_PORT", "5432"),
        "SUPABASE_DB_NAME": os.getenv("SUPABASE_DB_NAME", "postgres"),
        "SUPABASE_DB_USER": os.getenv("SUPABASE_DB_USER", "postgres"),
        "SUPABASE_DB_PASSWORD": os.getenv("SUPABASE_DB_PASSWORD", ""),
        
        # Database connection variables (DB_ prefix)
        "DB_HOST": os.getenv("DB_HOST", ""),
        "DB_PORT": os.getenv("DB_PORT", "5432"),
        "DB_NAME": os.getenv("DB_NAME", "postgres"),
        "DB_USER": os.getenv("DB_USER", "postgres"),
        "DB_PASSWORD": os.getenv("DB_PASSWORD", ""),
    }
    
    # Cross-fill missing variables
    if not env_vars["DB_HOST"] and env_vars["SUPABASE_DB_HOST"]:
        env_vars["DB_HOST"] = env_vars["SUPABASE_DB_HOST"]
    
    if not env_vars["SUPABASE_DB_HOST"] and env_vars["DB_HOST"]:
        env_vars["SUPABASE_DB_HOST"] = env_vars["DB_HOST"]
    
    # Same for other pairs of variables...
    
    return env_vars
```

### 3.3 Enhanced Connection Pool Implementation

#### 3.3.1 Key Features to Implement

1. **Graceful DNS Fallback**: Try hostname first, fall back to IP if DNS fails
2. **Improved Error Handling**: More detailed error logging and recovery
3. **Environment Variable Consistency**: Support for both DB_* and SUPABASE_DB_* patterns
4. **Health Checking**: Periodic connection validation and pool reinitialization

#### 3.3.2 Reference: Existing ConnectionPoolWrapper

From [`connection_pool_wrapper.py`](line:1-225):
```python
# Key method from current implementation
def _initialize_pool(self) -> None:
    """Initialize the connection pool with configuration parameters."""
    try:
        if self.pool is not None:
            # Close existing pool if it exists
            self.pool.closeall()
            
        self.pool = psycopg2.pool.ThreadedConnectionPool(
            self.min_conn,
            self.max_conn,
            host=self.config.get('host'),
            dbname=self.config.get('dbname'),
            user=self.config.get('user'),
            password=self.config.get('password'),
            port=self.config.get('port'),
            connect_timeout=self.connection_timeout
        )
        self.pool_initialized = True
        self.active_connections = 0
        logger.info("Connection pool initialized with min=%d, max=%d connections", 
                   self.min_conn, self.max_conn)
    except Exception as e:
        self.pool_initialized = False
        logger.error("Failed to initialize connection pool: %s", str(e))
        raise
```

#### 3.3.3 Example: Enhanced Connection Pool Initialization

```python
def _initialize_pool(self) -> None:
    """Initialize connection pool with hostname/IP fallback."""
    try:
        if self.pool is not None:
            # Close existing pool if it exists
            self.pool.closeall()
            
        # Create a copy of the config for connection
        conn_params = {
            'host': self.host,
            'dbname': self.config.get('dbname'),
            'user': self.config.get('user'),
            'password': self.config.get('password'),
            'port': self.config.get('port'),
            'connect_timeout': self.connection_timeout
        }
        
        try:
            # First try connecting with hostname
            self.pool = psycopg2.pool.ThreadedConnectionPool(
                self.min_conn,
                self.max_conn,
                **conn_params
            )
            self.pool_initialized = True
            self.active_connections = 0
            logger.info("Connection pool initialized with hostname %s", self.host)
        except (psycopg2.OperationalError, socket.gaierror) as e:
            # If hostname connection fails and we have an IP address, try with IP
            if self.ip_address:
                logger.warning("Hostname connection failed (%s), trying with IP %s", 
                               str(e), self.ip_address)
                conn_params['host'] = self.ip_address
                self.pool = psycopg2.pool.ThreadedConnectionPool(
                    self.min_conn,
                    self.max_conn,
                    **conn_params
                )
                self.pool_initialized = True
                self.active_connections = 0
                logger.info("Connection pool initialized with IP %s", self.ip_address)
            else:
                # Re-raise the exception if we don't have an IP fallback
                logger.error("Failed to initialize connection pool: %s", str(e))
                raise
    except Exception as e:
        self.pool_initialized = False
        logger.error("Failed to initialize connection pool: %s", str(e))
        raise
```

### 3.4 Update Config.py Implementation

#### 3.4.1 Current get_db_config() Function

From [`config.py`](line:615-633):
```python
def get_db_config() -> Dict[str, Any]:
    """
    Get database configuration with pool settings.
    
    Returns:
        Dict containing database configuration parameters
    """
    return {
        'host': os.environ.get('DB_HOST', 'localhost'),
        'dbname': os.environ.get('DB_NAME', 'cryoprotect'),
        'user': os.environ.get('DB_USER', 'postgres'),
        'password': os.environ.get('DB_PASSWORD', 'postgres'),
        'port': int(os.environ.get('DB_PORT', 5432)),
        'min_connections': int(os.environ.get('DB_MIN_CONNECTIONS', 1)),
        'max_connections': int(os.environ.get('DB_MAX_CONNECTIONS', 10)),
        'health_check_interval': int(os.environ.get('DB_HEALTH_CHECK_INTERVAL', 60)),
        'connection_timeout': int(os.environ.get('DB_CONNECTION_TIMEOUT', 30)),
    }
```

#### 3.4.2 Updated get_db_config() Function

```python
def get_db_config() -> Dict[str, Any]:
    """
    Get database configuration with pool settings.
    Supports both hostname and IP address connection.
    
    Returns:
        Dict containing database configuration parameters
    """
    # Load environment variables to ensure we have the latest
    load_dotenv()
    
    # Get host and IP
    host = os.environ.get('DB_HOST', os.environ.get('SUPABASE_DB_HOST', 'localhost'))
    ip = os.environ.get('SUPABASE_DB_IP')
    
    # Always use host in config, but the EnhancedConnectionPoolWrapper will
    # handle fallback to IP if hostname resolution fails
    return {
        'host': host,
        'dbname': os.environ.get('DB_NAME', os.environ.get('SUPABASE_DB_NAME', 'postgres')),
        'user': os.environ.get('DB_USER', os.environ.get('SUPABASE_DB_USER', 'postgres')),
        'password': os.environ.get('DB_PASSWORD', os.environ.get('SUPABASE_DB_PASSWORD', 'postgres')),
        'port': int(os.environ.get('DB_PORT', os.environ.get('SUPABASE_DB_PORT', 5432))),
        'min_connections': int(os.environ.get('DB_MIN_CONNECTIONS', 1)),
        'max_connections': int(os.environ.get('DB_MAX_CONNECTIONS', 10)),
        'health_check_interval': int(os.environ.get('DB_HEALTH_CHECK_INTERVAL', 60)),
        'connection_timeout': int(os.environ.get('DB_CONNECTION_TIMEOUT', 30)),
        'ip_address': ip,  # Optional IP address for fallback
    }
```

### 3.5 Hosts File Update Implementation (Optional)

#### 3.5.1 Function to Update Hosts File

```python
def update_hosts_file(hostname: str, ip_address: str) -> bool:
    """Update the hosts file to map hostname to IP."""
    # Determine hosts file location based on OS
    if platform.system() == "Windows":
        hosts_path = r"C:\Windows\System32\drivers\etc\hosts"
    else:
        hosts_path = "/etc/hosts"
    
    try:
        # Read current hosts file
        with open(hosts_path, 'r') as f:
            hosts_content = f.readlines()
        
        # Check if hostname is already in hosts file
        hostname_exists = False
        for i, line in enumerate(hosts_content):
            if not line.strip() or line.strip().startswith('#'):
                continue
                
            if hostname in line:
                # Update the line
                hosts_content[i] = f"{ip_address}\t{hostname}\n"
                hostname_exists = True
                break
        
        # Add hostname if it doesn't exist
        if not hostname_exists:
            api_hostname = hostname.replace("db.", "")
            hosts_content.append(f"\n# Added by CryoProtect v2 Supabase Connection Fix\n")
            hosts_content.append(f"{ip_address}\t{hostname}\n")
            hosts_content.append(f"{ip_address}\t{api_hostname}\n")
        
        # Write updated hosts file to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
            temp_path = temp_file.name
            temp_file.writelines(hosts_content)
        
        # Copy the temporary file to the hosts file with admin privileges
        # OS-specific implementation...
        
        return True
    except Exception as e:
        logger.error(f"Failed to update hosts file: {str(e)}")
        return False
```

## 4. Testing and Verification

### 4.1 Test Components

1. **DNS Resolution Test**: Verify hostname resolves to IP
2. **Direct Connection Test**: Test direct psycopg2 connection
3. **Connection Pool Test**: Test the enhanced connection pool
4. **Integration Test**: Verify with the existing verify_imported_data.py script

### 4.2 Example Test Script

```python
#!/usr/bin/env python3
"""
Test script to validate Supabase database connection
"""

import os
import socket
import psycopg2
import psycopg2.extras
import logging
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_dns_resolution():
    """Test DNS resolution for the Supabase database host."""
    load_dotenv()
    
    # Get database host from environment variables
    db_host = os.getenv('SUPABASE_DB_HOST', os.getenv('DB_HOST', ''))
    
    if not db_host:
        logger.error("No database host found in environment variables")
        return False
    
    logger.info(f"Testing DNS resolution for {db_host}...")
    
    try:
        ip_address = socket.gethostbyname(db_host)
        logger.info(f"Successfully resolved {db_host} to {ip_address}")
        return True
    except socket.gaierror as e:
        logger.error(f"Failed to resolve {db_host}: {str(e)}")
        
        # Check if we have a direct IP in environment
        db_ip = os.getenv('SUPABASE_DB_IP', '')
        if db_ip:
            logger.info(f"Found direct IP in environment: {db_ip}")
            return True
        
        return False

def test_direct_connection():
    """Test direct PostgreSQL connection."""
    # Implementation...

def test_enhanced_connection_pool():
    """Test the enhanced connection pool wrapper."""
    # Implementation...

def main():
    """Run all connection tests."""
    print("CryoProtect v2 - Supabase Connection Test")
    
    # Test DNS resolution
    dns_ok = test_dns_resolution()
    
    # Test direct connection
    direct_ok = test_direct_connection()
    
    # Test enhanced connection pool
    pool_ok = test_enhanced_connection_pool()
    
    # Print summary
    print("\nTest Summary:")
    print(f"  DNS Resolution:       {'✓' if dns_ok else '✗'}")
    print(f"  Direct Connection:    {'✓' if direct_ok else '✗'}")
    print(f"  Enhanced Pool:        {'✓' if pool_ok else '✗'}")
    
    if dns_ok and direct_ok and pool_ok:
        print("\nSUCCESS: All tests passed!")
        return 0
    else:
        print("\nWARNING: Some tests failed.")
        return 1

if __name__ == "__main__":
    main()
```

## 5. Implementation Task Breakdown

### 5.1 Task List

1. **Create IP Resolution Module**
   - Implement multiple resolution methods
   - Add fallback mechanisms
   - Handle error cases gracefully

2. **Update Environment Configuration**
   - Modify .env loading
   - Ensure consistent variable patterns
   - Add support for IP-based connection

3. **Enhance Connection Pool**
   - Create enhanced connection pool wrapper
   - Add DNS failure recovery
   - Implement graceful degradation

4. **Update Config Module**
   - Enhance get_db_config() function
   - Add support for both variable patterns
   - Include IP address support

5. **Create Testing Tools**
   - Implement validation script
   - Add diagnostic functions
   - Create easy-to-run tests

### 5.2 Implementation Order

1. First implement the IP resolution module (highest priority)
2. Then update environment configuration to include IP address
3. Next enhance the connection pool to handle DNS failures
4. After that update config.py to use the new environment variables
5. Finally create testing tools to validate the solution

## 6. Additional Resources

### 6.1 Code References

- **Current Connection Pool**: [connection_pool_wrapper.py](lines 1-225)
- **Database Config**: [config.py](lines 615-633)
- **Verification Script**: [verify_imported_data.py](lines 1-481)
- **DNS Fix Attempt**: [fix_dns_mcp.py](lines 1-221)
- **Direct Connection Attempt**: [direct_postgres_connect.py](lines 1-570)

### 6.2 Related Documentation

- **Supabase API**: [https://supabase.io/docs/reference/javascript/initializing](https://supabase.io/docs/reference/javascript/initializing)
- **Psycopg2 Connection**: [https://www.psycopg.org/docs/module.html#psycopg2.connect](https://www.psycopg.org/docs/module.html#psycopg2.connect)
- **Connection Pooling**: [https://www.psycopg.org/docs/pool.html](https://www.psycopg.org/docs/pool.html)

## 7. Project Integration Notes

### 7.1 Backward Compatibility Requirements

- The solution must maintain compatibility with existing code
- It should work without requiring significant changes to application code
- The original connection_pool_wrapper.py should still function (but new version preferred)

### 7.2 Deployment Considerations

- The solution should work on both Windows and Unix-like systems
- No additional dependencies should be introduced if possible
- Changes should be documented for future reference

### 7.3 Error Handling Guidelines

- All errors should be logged with sufficient detail
- User-friendly error messages should be provided
- Fallback mechanisms should be implemented where possible

## 8. Implementation Example

Here's a minimal implementation example for the main connection fix:

```python
#!/usr/bin/env python3
"""
Enhanced database connection module for CryoProtect v2
"""

import os
import socket
import logging
from typing import Dict, Any, Optional
import psycopg2
import psycopg2.pool
import psycopg2.extras
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class EnhancedDatabaseConnection:
    """Enhanced database connection with DNS fallback."""
    
    def __init__(self):
        # Get connection parameters
        self.host = os.getenv('SUPABASE_DB_HOST', os.getenv('DB_HOST', 'localhost'))
        self.port = int(os.getenv('SUPABASE_DB_PORT', os.getenv('DB_PORT', '5432')))
        self.dbname = os.getenv('SUPABASE_DB_NAME', os.getenv('DB_NAME', 'postgres'))
        self.user = os.getenv('SUPABASE_DB_USER', os.getenv('DB_USER', 'postgres'))
        self.password = os.getenv('SUPABASE_DB_PASSWORD', os.getenv('DB_PASSWORD', ''))
        
        # Try to resolve IP address
        self.ip_address = os.getenv('SUPABASE_DB_IP', '')
        if not self.ip_address and self.host:
            try:
                self.ip_address = socket.gethostbyname(self.host)
                logger.info(f"Resolved {self.host} to {self.ip_address}")
            except socket.gaierror:
                logger.warning(f"Could not resolve {self.host}")
    
    def get_connection(self):
        """Get a database connection with fallback to IP."""
        try:
            # Try connecting with hostname
            conn = psycopg2.connect(
                host=self.host,
                port=self.port,
                dbname=self.dbname,
                user=self.user,
                password=self.password,
                connect_timeout=10
            )
            logger.info(f"Connected to {self.host} successfully")
            return conn
        except (psycopg2.OperationalError, socket.gaierror) as e:
            logger.warning(f"Connection to {self.host} failed: {str(e)}")
            
            # Try with IP if available
            if self.ip_address:
                logger.info(f"Trying connection with IP {self.ip_address}")
                try:
                    conn = psycopg2.connect(
                        host=self.ip_address,
                        port=self.port,
                        dbname=self.dbname,
                        user=self.user,
                        password=self.password,
                        connect_timeout=10
                    )
                    logger.info(f"Connected to {self.ip_address} successfully")
                    return conn
                except Exception as e2:
                    logger.error(f"Connection to IP {self.ip_address} failed: {str(e2)}")
            
            # Re-raise the original exception
            raise e

# Example usage
if __name__ == "__main__":
    db = EnhancedDatabaseConnection()
    try:
        conn = db.get_connection()
        with conn.cursor() as cursor:
            cursor.execute("SELECT 1")
            result = cursor.fetchone()
            print(f"Connection test result: {result[0]}")
        conn.close()
        print("Connection test passed!")
    except Exception as e:
        print(f"Connection test failed: {str(e)}")
```

## 9. Conclusion

Implementing this comprehensive solution will resolve the database connection issues in the CryoProtect v2 project. The multi-faceted approach addresses DNS resolution, configuration inconsistency, and connection pooling issues, providing a robust and reliable solution.

By following this guide, Roo code agents can implement a solution that:

1. Resolves IP addresses using multiple methods
2. Ensures consistent environment variables
3. Enhances the connection pool to handle failures gracefully
4. Includes comprehensive testing and validation

This will ensure reliable database connections for the application, enabling successful database population and verification.