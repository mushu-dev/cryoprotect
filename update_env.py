#!/usr/bin/env python3
from dotenv import load_dotenv, set_key

# Load current .env file
load_dotenv()

# Update with correct values
set_key('.env', 'SUPABASE_URL', 'https://tsdlmynydfuypiugmkev.supabase.co')
set_key('.env', 'SUPABASE_DB_HOST', 'db.tsdlmynydfuypiugmkev.supabase.co')

print("Updated .env file with correct hostnames")
