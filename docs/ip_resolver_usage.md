# IP Resolution Module Documentation

## Overview

The IP Resolution module (`ip_resolver.py`) provides robust hostname-to-IP resolution capabilities for the CryoProtect v2 project. It was specifically designed to address DNS resolution issues with Supabase database connections, but can be used for any hostname resolution needs.

## Features

- **Multiple Resolution Methods**: Uses a fallback chain of resolution methods:
  1. Standard DNS resolution via `socket.gethostbyname()`
  2. Alternative DNS servers (8.8.8.8, 1.1.1.1) using OS-specific commands
  3. Placeholder for MCP-based resolution (for future implementation)
  4. Placeholder for heuristic-based resolution (for future implementation)

- **OS-Specific Implementation**: Handles differences between Windows and Unix systems:
  - Windows: Uses `nslookup` command
  - Unix: Uses `dig` command

- **Environment Variable Management**: Can automatically update `.env` file with resolved IP addresses

- **Comprehensive Logging**: Detailed logging of all resolution attempts and results

- **Robust Error Handling**: Gracefully handles all error conditions

## Usage

### Basic Usage

```python
from ip_resolver import resolve_ip_address

# Resolve a hostname to an IP address
ip = resolve_ip_address("db.tsdlmynydfuypiugmkev.supabase.co")
if ip:
    print(f"Resolved to: {ip}")
else:
    print("Resolution failed")
```

### Updating Environment Variables

```python
from ip_resolver import resolve_and_update_env

# Resolve hostname and update .env file
result = resolve_and_update_env("db.tsdlmynydfuypiugmkev.supabase.co")
if result['ip']:
    print(f"Resolved to: {result['ip']}")
    if result['updated']:
        print("Updated .env file")
else:
    print("Resolution failed")
```

### Command Line Usage

The module can also be run directly from the command line:

```bash
python ip_resolver.py [hostname]
```

If no hostname is provided, it defaults to "db.tsdlmynydfuypiugmkev.supabase.co".

## Testing

Two test files are provided:

1. **Simple Test Script**: `test_ip_resolver.py` - A simple script to test the module's functionality
   ```bash
   python test_ip_resolver.py [hostname]
   ```

2. **Unit Tests**: `tests/unit/test_ip_resolver.py` - Comprehensive unit tests for the module
   ```bash
   python -m unittest tests/unit/test_ip_resolver.py
   ```

## Integration with Connection Pool

To use the IP resolver with the connection pool wrapper:

1. First resolve the hostname to an IP address:
   ```python
   from ip_resolver import resolve_ip_address
   
   hostname = "db.tsdlmynydfuypiugmkev.supabase.co"
   ip_address = resolve_ip_address(hostname)
   ```

2. Update your database configuration to include the IP address:
   ```python
   db_config = {
       'host': hostname,  # Keep the hostname for standard connection attempts
       'ip_address': ip_address,  # Add the IP address for fallback
       # Other database configuration parameters...
   }
   ```

3. Use the enhanced connection pool wrapper that supports IP fallback.

## Error Handling

The module includes comprehensive error handling for various failure scenarios:

- Socket errors during standard DNS resolution
- Command execution errors during alternative DNS resolution
- Timeout errors during command execution
- File I/O errors during .env file updates

All errors are logged with appropriate severity levels and meaningful messages.

## Future Enhancements

The module includes placeholder implementations for two additional resolution methods:

1. **MCP Resolution**: Using the Supabase MCP to get project information including IP addresses
2. **Heuristic Resolution**: Using patterns in IP assignment as a fallback method

These methods can be implemented in the future to further enhance the module's capabilities.