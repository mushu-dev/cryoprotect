#!/usr/bin/env python3
"""
This script identifies users without a team assignment and assigns them to a
default team. It's designed to ensure all users have proper team association
for permissions and resource access.
"""

import os
import sys
import logging
import psycopg2
import psycopg2.extras
from dotenv import load_dotenv
from datetime import datetime
import uuid
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("team_assignment.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST")
DB_PORT = os.getenv("SUPABASE_DB_PORT")
DB_NAME = os.getenv("SUPABASE_DB_NAME")
DB_USER = os.getenv("SUPABASE_DB_USER")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD")

def get_db_connection():
    """Get a direct database connection using psycopg2."""
    try:
        conn = psycopg2.connect(
            host=DB_HOST,
            port=DB_PORT,
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD
        )
        logger.info("Connected to database using direct PostgreSQL connection")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def get_or_create_default_team():
    """
    Get the default team or create it if it doesn't exist.
    
    Returns:
        dict: Team information with id and name
    """
    conn = get_db_connection()
    if not conn:
        return None
        
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    team = None
    
    try:
        # Look for the default team
        cursor.execute("""
            SELECT id, name, description
            FROM teams
            WHERE name = 'Default Team'
        """)
        
        team = cursor.fetchone()
        
        if team:
            logger.info(f"Found default team: {team['name']} (ID: {team['id']})")
        else:
            # Create the default team
            team_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO teams (id, name, description, created_at, updated_at)
                VALUES (%s, %s, %s, NOW(), NOW())
                RETURNING id, name, description
            """, (
                team_id,
                'Default Team',
                'Default team for users without team assignment'
            ))
            
            team = cursor.fetchone()
            conn.commit()
            logger.info(f"Created default team: {team['name']} (ID: {team['id']})")
            
    except Exception as e:
        conn.rollback()
        logger.error(f"Error getting or creating default team: {e}")
    
    cursor.close()
    conn.close()
    return team

def get_users_without_team():
    """
    Get users who are not assigned to any team.
    
    Returns:
        list: User records without team assignment
    """
    conn = get_db_connection()
    if not conn:
        return []
        
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    users = []
    
    try:
        cursor.execute("""
            SELECT u.id, u.email, p.display_name
            FROM auth.users u
            LEFT JOIN user_profile p ON u.id = p.auth_user_id
            WHERE u.id NOT IN (
                SELECT user_id FROM team_members
            )
            AND u.confirmed_at IS NOT NULL
        """)
        
        users = cursor.fetchall()
        logger.info(f"Found {len(users)} users without team assignment")
        
    except Exception as e:
        logger.error(f"Error getting users without team: {e}")
    
    cursor.close()
    conn.close()
    return users

def assign_user_to_team(user_id, team_id, role="member", dry_run=False):
    """
    Assign a user to a team.
    
    Args:
        user_id: User ID to assign
        team_id: Team ID to assign to
        role: Role in the team (default: member)
        dry_run: If True, don't make actual changes
        
    Returns:
        bool: True if successful, False otherwise
    """
    conn = get_db_connection()
    if not conn:
        return False
        
    cursor = conn.cursor()
    success = False
    
    try:
        if dry_run:
            logger.info(f"DRY RUN: Would assign user {user_id} to team {team_id} with role '{role}'")
            success = True
        else:
            cursor.execute("""
                INSERT INTO team_members (team_id, user_id, role, created_at, updated_at)
                VALUES (%s, %s, %s, NOW(), NOW())
                ON CONFLICT (team_id, user_id) DO NOTHING
            """, (team_id, user_id, role))
            
            if cursor.rowcount > 0:
                logger.info(f"Assigned user {user_id} to team {team_id} with role '{role}'")
            else:
                logger.warning(f"User {user_id} already assigned to team {team_id}")
                
            conn.commit()
            success = True
            
    except Exception as e:
        conn.rollback()
        logger.error(f"Error assigning user {user_id} to team {team_id}: {e}")
    
    cursor.close()
    conn.close()
    return success

def main():
    parser = argparse.ArgumentParser(description="Assign users without a team to the default team.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing them")
    parser.add_argument("--role", default="member", choices=["member", "admin", "owner"], 
                        help="Role to assign to users (default: member)")
    args = parser.parse_args()
    
    logger.info(f"Starting team assignment {'in dry-run mode' if args.dry_run else 'in execution mode'}")
    
    # Get or create the default team
    default_team = get_or_create_default_team()
    if not default_team:
        logger.error("Failed to get or create default team")
        return False
        
    # Get users without a team
    users = get_users_without_team()
    if not users:
        logger.info("No users found without team assignment")
        return True
        
    # Assign users to the default team
    assigned_count = 0
    for user in users:
        user_id = user['id']
        
        if args.dry_run:
            logger.info(f"DRY RUN: Would assign user {user.get('email', user_id)} to '{default_team['name']}'")
            assigned_count += 1
        else:
            if assign_user_to_team(user_id, default_team['id'], args.role):
                assigned_count += 1
    
    # Print summary
    print("\n" + "=" * 60)
    print("Team Assignment Summary")
    print("=" * 60)
    print(f"Default team: {default_team['name']} (ID: {default_team['id']})")
    print(f"Total users processed: {len(users)}")
    print(f"Users assigned: {assigned_count}")
    print("=" * 60)
    
    logger.info("Team assignment complete")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)