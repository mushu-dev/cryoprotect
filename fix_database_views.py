#!/usr/bin/env python3
"""
fix_database_views.py - Create missing database views with plural names

This script creates the missing database views with plural names to match the updated API code.
The API code has been updated to use plural table names, but the database views haven't been
created or renamed yet.

Key actions:
1. Create molecules_with_properties view (previously molecule_with_properties)
2. Create mixtures_with_components view (previously mixture_with_components)
3. Verify that the views have been created correctly

Usage:
    python fix_database_views.py
"""

import os
import sys
import logging
from datetime import datetime
from pathlib import Path
import json

# Configure logging
logs_dir = Path("logs")
logs_dir.mkdir(exist_ok=True)

log_file = logs_dir / f"database_views_fix_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Import Supabase client utility
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from api.utils import get_supabase_client, execute_with_retry
except ImportError:
    logger.error("Failed to import Supabase client utility. Make sure you're running this script from the project root.")
    sys.exit(1)

def check_view_exists(view_name):
    """Check if a view exists in the database."""
    try:
        supabase = get_supabase_client()
        
        # Query the information_schema.views table to check if the view exists
        response = supabase.from_("information_schema.views").select("*").eq("table_name", view_name).execute()
        
        return len(response.data) > 0
    except Exception as e:
