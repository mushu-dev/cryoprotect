# CryoProtect v2 - Supabase Direct Connection Workaround

This directory contains workarounds for DNS resolution issues with Supabase.

## Issue

The issue is that the database hostname (`db.tsdlmynydfuypiugmkev.supabase.co`) cannot be resolved
through DNS. This causes connections to fail when trying to connect to the database directly.

## Solution

The workaround is to use direct IP connections instead of hostnames:

1. Use the direct IP address for the API hostname
2. Add a `Host` header to requests to ensure SSL certificate validation works
3. Create a helper class for making requests to Supabase using direct IP

## Files

- `direct_connection_fix.py`: Script that sets up the workaround
- `supabase_direct_connection.py`: Helper class for connecting to Supabase
- `test_direct_connection.py`: Script for testing the workaround

## Usage

To use the workaround in your code:

```python
from supabase_direct_connection import SupabaseDirectConnection

# Create connection
supabase = SupabaseDirectConnection()

# Get data
molecules = supabase.get_data("molecules", {"limit": 10})
print(f"Retrieved {len(molecules)} molecules")
```

## Important Notes

- This is a temporary workaround until DNS resolution issues are fixed
- The direct IP address may change in the future
- You should periodically try to connect using the original hostname to see if the issue is resolved
