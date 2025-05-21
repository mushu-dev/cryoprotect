try:
    import jwt
    print(f"JWT module is available: {jwt}")
    print(f"JWT version: {jwt.__version__ if hasattr(jwt, '__version__') else 'unknown'}")
except ImportError as e:
    print(f"Failed to import jwt: {e}")
