"""
CryoProtect v2 - Service Role Helper

This module provides helper functions for connecting to Supabase using the service role
and managing user authentication and profiles.
"""

import os
import uuid
import logging
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase connection
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SUPABASE_USER = os.getenv("SUPABASE_USER")
SUPABASE_PASSWORD = os.getenv("SUPABASE_PASSWORD")

# Cache for user ID and profile ID
_user_id_cache = None
_profile_id_cache = None

def get_supabase_client() -> Client:
    """
    Get a Supabase client using the service role key.
    
    Returns:
        Client: Supabase client with service role permissions
    """
    if not SUPABASE_URL or not SUPABASE_KEY:
        raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
    
    return create_client(SUPABASE_URL, SUPABASE_KEY)

def get_user_id() -> str:
    """
    Get the user ID of the authenticated user.
    
    Returns:
        str: User ID or None if not authenticated
    """
    global _user_id_cache
    
    # Return cached user ID if available
    if _user_id_cache:
        return _user_id_cache
    
    if not SUPABASE_USER or not SUPABASE_PASSWORD:
        logger.warning("SUPABASE_USER and SUPABASE_PASSWORD must be set for authentication")
        return None
    
    try:
        supabase = get_supabase_client()
        
        # Sign in with email and password
        response = supabase.auth.sign_in_with_password({
            "email": SUPABASE_USER,
            "password": SUPABASE_PASSWORD
        })
        
        if hasattr(response, 'user') and response.user:
            _user_id_cache = response.user.id
            logger.info(f"Authenticated as {SUPABASE_USER}")
            return _user_id_cache
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Authentication error: {response.error}")
        
        return None
    except Exception as e:
        logger.error(f"Error getting user ID: {str(e)}")
        return None

def ensure_user_profile(supabase: Client = None) -> str:
    """
    Ensure a user profile exists for the authenticated user.
    
    Args:
        supabase (Client, optional): Supabase client. If None, a new client will be created.
        
    Returns:
        str: Profile ID or None if profile creation failed
    """
    global _profile_id_cache
    
    # Return cached profile ID if available
    if _profile_id_cache:
        return _profile_id_cache
    
    # Get user ID
    user_id = get_user_id()
    if not user_id:
        logger.warning("No authenticated user ID available")
        return None
    
    # Create Supabase client if not provided
    if not supabase:
        supabase = get_supabase_client()
    
    try:
        # Check if profile already exists
        response = supabase.table("user_profile").select("*").eq("auth_user_id", user_id).execute()
        
        if hasattr(response, 'data') and response.data:
            _profile_id_cache = response.data[0]["id"]
            logger.info(f"User profile already exists with ID: {_profile_id_cache}")
            return _profile_id_cache
        
        # Create new profile
        profile_id = str(uuid.uuid4())
        profile_data = {
            "id": profile_id,
            "auth_user_id": user_id,
            "display_name": SUPABASE_USER.split('@')[0] if SUPABASE_USER else "CryoProtect User",
            "email": SUPABASE_USER,
            "affiliation": "CryoProtect Project",
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat()
        }
        
        response = supabase.table("user_profile").insert(profile_data).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error creating user profile: {response.error}")
            return None
        
        _profile_id_cache = profile_id
        logger.info(f"Created user profile with ID: {profile_id}")
        return profile_id
    
    except Exception as e:
        logger.error(f"Error ensuring user profile: {str(e)}")
        return None

def get_project_id(supabase: Client = None, project_id: str = None) -> str:
    """
    Get a project ID, either from parameter or by finding the first available project.
    
    Args:
        supabase (Client, optional): Supabase client. If None, a new client will be created.
        project_id (str, optional): Specific project ID to use.
        
    Returns:
        str: Project ID or None if no project found
    """
    # Create Supabase client if not provided
    if not supabase:
        supabase = get_supabase_client()
    
    if project_id:
        # Verify the project exists
        response = supabase.table("projects").select("id").eq("id", project_id).execute()
        if hasattr(response, 'data') and response.data:
            logger.info(f"Using specified project ID: {project_id}")
            return project_id
        else:
            logger.warning(f"Specified project ID {project_id} not found. Will try to find another project.")
    
    # Try to find any project
    response = supabase.table("projects").select("id").execute()
    if hasattr(response, 'data') and response.data:
        project_id = response.data[0]["id"]
        logger.info(f"Found existing project ID: {project_id}")
        return project_id
    
    logger.warning("No projects found. Data will be created without project association.")
    return None

if __name__ == "__main__":
    # Test the helper functions
    supabase = get_supabase_client()
    user_id = get_user_id()
    profile_id = ensure_user_profile(supabase)
    project_id = get_project_id(supabase)
    
    print(f"User ID: {user_id}")
    print(f"Profile ID: {profile_id}")
    print(f"Project ID: {project_id}")
