#!/usr/bin/env python3
"""
Script to implement proper team-based access control with RLS.
This script applies the SQL migration and tests the functionality.
"""

import os
import sys
import json
import logging
import uuid
from datetime import datetime
import argparse
from psycopg2.extras import RealDictCursor
import db_utils

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_migration():
    """Execute the migration to implement team-based access control."""
    sql_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        'migrations', 
        '030_team_based_access_control.sql'
    )
    
    try:
        if db_utils.execute_sql_file(sql_file):
            logger.info("Successfully executed migration to implement team-based access control")
            return True
        else:
            logger.error("Failed to execute migration")
            return False
    except Exception as e:
        logger.error(f"Error executing migration: {e}")
        return False

def check_teams_structure():
    """Check if the teams and team_members tables exist with the correct structure."""
    try:
        # Check if teams table exists
        teams_query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'teams'
        );
        """
        teams_exists = db_utils.execute_query(teams_query)[0][0]
        
        # Check if team_members table exists
        members_query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'team_members'
        );
        """
        members_exists = db_utils.execute_query(members_query)[0][0]
        
        # Check if team_invitations table exists
        invitations_query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'team_invitations'
        );
        """
        invitations_exists = db_utils.execute_query(invitations_query)[0][0]
        
        # Check team columns if exists
        team_columns = []
        if teams_exists:
            columns_query = """
            SELECT column_name, data_type
            FROM information_schema.columns
            WHERE table_schema = 'public'
            AND table_name = 'teams';
            """
            team_columns = db_utils.execute_query(columns_query, cursor_factory=RealDictCursor)
        
        # Check team_members columns if exists
        member_columns = []
        if members_exists:
            columns_query = """
            SELECT column_name, data_type
            FROM information_schema.columns
            WHERE table_schema = 'public'
            AND table_name = 'team_members';
            """
            member_columns = db_utils.execute_query(columns_query, cursor_factory=RealDictCursor)
        
        return {
            'teams_exists': teams_exists,
            'team_members_exists': members_exists,
            'team_invitations_exists': invitations_exists,
            'team_columns': team_columns,
            'team_member_columns': member_columns
        }
    except Exception as e:
        logger.error(f"Error checking teams structure: {e}")
        return {
            'teams_exists': False,
            'team_members_exists': False,
            'team_invitations_exists': False,
            'error': str(e)
        }

def create_team_tables():
    """Create the team-related tables if they don't exist."""
    try:
        # Check if teams table exists
        structure = check_teams_structure()
        
        if not structure['teams_exists']:
            # Create teams table
            teams_query = """
            CREATE TABLE public.teams (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                name TEXT NOT NULL,
                description TEXT,
                created_by UUID NOT NULL REFERENCES auth.users(id),
                created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
            );
            """
            db_utils.execute_query(teams_query, fetch=False)
            logger.info("Created teams table")
        else:
            logger.info("Teams table already exists")
        
        if not structure['team_members_exists']:
            # Create team_members table
            members_query = """
            CREATE TABLE public.team_members (
                team_id UUID NOT NULL REFERENCES public.teams(id) ON DELETE CASCADE,
                user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
                role TEXT NOT NULL DEFAULT 'read',
                joined_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                PRIMARY KEY (team_id, user_id),
                CONSTRAINT valid_role CHECK (role IN ('admin', 'write', 'read'))
            );
            """
            db_utils.execute_query(members_query, fetch=False)
            logger.info("Created team_members table")
        else:
            logger.info("Team_members table already exists")
        
        # Add team_id to molecules and mixtures tables if needed
        add_team_id_query = """
        DO $$
        BEGIN
            -- Check if team_id column exists in molecules table
            IF NOT EXISTS (
                SELECT FROM information_schema.columns 
                WHERE table_schema = 'public' 
                AND table_name = 'molecules' 
                AND column_name = 'team_id'
            ) THEN
                ALTER TABLE public.molecules ADD COLUMN team_id UUID REFERENCES public.teams(id);
                CREATE INDEX idx_molecules_team_id ON public.molecules(team_id);
                RAISE NOTICE 'Added team_id column to molecules table';
            END IF;
            
            -- Check if team_id column exists in mixtures table
            IF NOT EXISTS (
                SELECT FROM information_schema.columns 
                WHERE table_schema = 'public' 
                AND table_name = 'mixtures' 
                AND column_name = 'team_id'
            ) THEN
                ALTER TABLE public.mixtures ADD COLUMN team_id UUID REFERENCES public.teams(id);
                CREATE INDEX idx_mixtures_team_id ON public.mixtures(team_id);
                RAISE NOTICE 'Added team_id column to mixtures table';
            END IF;
        END $$;
        """
        db_utils.execute_query(add_team_id_query, fetch=False)
        
        return True
    except Exception as e:
        logger.error(f"Error creating team tables: {e}")
        return False

def test_team_functions():
    """Test the team-related functions."""
    try:
        # Test auth.is_team_member function
        is_team_member_query = """
        SELECT pg_get_functiondef('auth.is_team_member(uuid)'::regprocedure);
        """
        is_team_member_result = db_utils.execute_query(is_team_member_query)[0][0]
        
        # Test auth.is_team_admin function
        is_team_admin_query = """
        SELECT pg_get_functiondef('auth.is_team_admin(uuid)'::regprocedure);
        """
        is_team_admin_result = db_utils.execute_query(is_team_admin_query)[0][0]
        
        # Test auth.get_user_teams function
        get_user_teams_query = """
        SELECT pg_get_functiondef('auth.get_user_teams()'::regprocedure);
        """
        get_user_teams_result = db_utils.execute_query(get_user_teams_query)[0][0]
        
        # Test auth.create_team function
        create_team_query = """
        SELECT pg_get_functiondef('auth.create_team(text,text)'::regprocedure);
        """
        create_team_result = db_utils.execute_query(create_team_query)[0][0]
        
        # Test auth.invite_to_team function
        invite_to_team_query = """
        SELECT pg_get_functiondef('auth.invite_to_team(uuid,text,text)'::regprocedure);
        """
        invite_to_team_result = db_utils.execute_query(invite_to_team_query)[0][0]
        
        # Test auth.accept_team_invitation function
        accept_invitation_query = """
        SELECT pg_get_functiondef('auth.accept_team_invitation(text)'::regprocedure);
        """
        accept_invitation_result = db_utils.execute_query(accept_invitation_query)[0][0]
        
        # Test auth.transfer_resource_to_team function
        transfer_resource_query = """
        SELECT pg_get_functiondef('auth.transfer_resource_to_team(text,uuid,uuid)'::regprocedure);
        """
        transfer_resource_result = db_utils.execute_query(transfer_resource_query)[0][0]
        
        return {
            'is_team_member': is_team_member_result is not None,
            'is_team_admin': is_team_admin_result is not None,
            'get_user_teams': get_user_teams_result is not None,
            'create_team': create_team_result is not None,
            'invite_to_team': invite_to_team_result is not None,
            'accept_team_invitation': accept_invitation_result is not None,
            'transfer_resource_to_team': transfer_resource_result is not None
        }
    except Exception as e:
        logger.error(f"Error testing team functions: {e}")
        return {
            'error': str(e)
        }

def test_rls_policies():
    """Test the RLS policies for team-based access control."""
    try:
        # Check RLS policies on teams table
        teams_policies_query = """
        SELECT policyname, cmd, qual
        FROM pg_policies
        WHERE schemaname = 'public'
        AND tablename = 'teams';
        """
        teams_policies = db_utils.execute_query(teams_policies_query, cursor_factory=RealDictCursor)
        
        # Check RLS policies on team_members table
        members_policies_query = """
        SELECT policyname, cmd, qual
        FROM pg_policies
        WHERE schemaname = 'public'
        AND tablename = 'team_members';
        """
        members_policies = db_utils.execute_query(members_policies_query, cursor_factory=RealDictCursor)
        
        # Check team-specific RLS policies on molecules table
        molecules_policies_query = """
        SELECT policyname, cmd, qual
        FROM pg_policies
        WHERE schemaname = 'public'
        AND tablename = 'molecules'
        AND policyname LIKE '%team%';
        """
        molecules_policies = db_utils.execute_query(molecules_policies_query, cursor_factory=RealDictCursor)
        
        # Check team-specific RLS policies on mixtures table
        mixtures_policies_query = """
        SELECT policyname, cmd, qual
        FROM pg_policies
        WHERE schemaname = 'public'
        AND tablename = 'mixtures'
        AND policyname LIKE '%team%';
        """
        mixtures_policies = db_utils.execute_query(mixtures_policies_query, cursor_factory=RealDictCursor)
        
        return {
            'teams_policies': teams_policies,
            'team_members_policies': members_policies,
            'molecules_team_policies': molecules_policies,
            'mixtures_team_policies': mixtures_policies
        }
    except Exception as e:
        logger.error(f"Error testing RLS policies: {e}")
        return {
            'error': str(e)
        }

def check_team_view():
    """Check if the team resources view exists and has the correct structure."""
    try:
        view_exists_query = """
        SELECT EXISTS (
            SELECT FROM information_schema.views 
            WHERE table_schema = 'auth' 
            AND table_name = 'team_resources'
        );
        """
        view_exists = db_utils.execute_query(view_exists_query)[0][0]
        
        view_columns = []
        if view_exists:
            view_columns_query = """
            SELECT column_name, data_type
            FROM information_schema.columns
            WHERE table_schema = 'auth'
            AND table_name = 'team_resources';
            """
            view_columns = db_utils.execute_query(view_columns_query, cursor_factory=RealDictCursor)
        
        return {
            'view_exists': view_exists,
            'view_columns': view_columns
        }
    except Exception as e:
        logger.error(f"Error checking team resources view: {e}")
        return {
            'view_exists': False,
            'error': str(e)
        }

def save_report(team_structure, function_tests, policy_tests, view_check, 
                filename_prefix="team_based_access_control_report"):
    """Save the report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    # Calculate success metrics
    structure_success = team_structure['teams_exists'] and team_structure['team_members_exists']
    function_success = all(v is True for k, v in function_tests.items() if k != 'error')
    policy_success = (
        len(policy_tests.get('teams_policies', [])) > 0 and
        len(policy_tests.get('team_members_policies', [])) > 0 and
        len(policy_tests.get('molecules_team_policies', [])) > 0 and
        len(policy_tests.get('mixtures_team_policies', [])) > 0
    )
    view_success = view_check.get('view_exists', False)
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'team_structure': team_structure,
        'function_tests': function_tests,
        'policy_tests': policy_tests,
        'view_check': view_check,
        'summary': {
            'structure_success': structure_success,
            'function_success': function_success,
            'policy_success': policy_success,
            'view_success': view_success,
            'overall_success': structure_success and function_success and policy_success and view_success
        }
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def create_test_data():
    """Create test data to validate team-based access control."""
    try:
        # Create a test team
        team_id = str(uuid.uuid4())
        team_query = f"""
        INSERT INTO public.teams (id, name, description, created_by)
        VALUES ('{team_id}', 'Test Team', 'Test team for validation', auth.uid())
        ON CONFLICT DO NOTHING
        RETURNING id;
        """
        team_result = db_utils.execute_query(team_query)
        if team_result:
            team_id = team_result[0][0]
            logger.info(f"Created test team with ID: {team_id}")
        else:
            logger.info("Test team already exists")
        
        # Create a test team member
        member_query = f"""
        INSERT INTO public.team_members (team_id, user_id, role)
        VALUES ('{team_id}', auth.uid(), 'admin')
        ON CONFLICT DO NOTHING;
        """
        db_utils.execute_query(member_query, fetch=False)
        logger.info("Added current user as team admin")
        
        # Create a test molecule assigned to the team
        molecule_id = str(uuid.uuid4())
        molecule_query = f"""
        INSERT INTO public.molecules (id, name, smiles, formula, team_id, created_by, is_public)
        VALUES (
            '{molecule_id}', 
            'Test Team Molecule', 
            'CC(=O)OC1=CC=CC=C1C(=O)O', 
            'C9H8O4', 
            '{team_id}', 
            auth.uid(),
            FALSE
        )
        ON CONFLICT DO NOTHING
        RETURNING id;
        """
        molecule_result = db_utils.execute_query(molecule_query)
        if molecule_result:
            molecule_id = molecule_result[0][0]
            logger.info(f"Created test molecule with ID: {molecule_id}")
        else:
            logger.info("Test molecule already exists")
        
        # Create a test mixture assigned to the team
        mixture_id = str(uuid.uuid4())
        mixture_query = f"""
        INSERT INTO public.mixtures (id, name, description, team_id, created_by, is_public)
        VALUES (
            '{mixture_id}', 
            'Test Team Mixture', 
            'Test mixture for team access control', 
            '{team_id}', 
            auth.uid(),
            FALSE
        )
        ON CONFLICT DO NOTHING
        RETURNING id;
        """
        mixture_result = db_utils.execute_query(mixture_query)
        if mixture_result:
            mixture_id = mixture_result[0][0]
            logger.info(f"Created test mixture with ID: {mixture_id}")
        else:
            logger.info("Test mixture already exists")
        
        # Create a test invitation
        invitation_query = f"""
        INSERT INTO public.team_invitations (team_id, email, role, invited_by)
        VALUES (
            '{team_id}', 
            'test@example.com', 
            'read', 
            auth.uid()
        )
        ON CONFLICT DO NOTHING
        RETURNING id, token;
        """
        invitation_result = db_utils.execute_query(invitation_query, cursor_factory=RealDictCursor)
        invitation_id = None
        invitation_token = None
        if invitation_result:
            invitation_id = invitation_result[0]['id']
            invitation_token = invitation_result[0]['token']
            logger.info(f"Created test invitation with ID: {invitation_id}")
        else:
            logger.info("Test invitation already exists")
        
        # Query the team resources view
        resources_query = f"""
        SELECT * FROM auth.team_resources
        WHERE team_id = '{team_id}';
        """
        resources = db_utils.execute_query(resources_query, cursor_factory=RealDictCursor)
        logger.info(f"Found {len(resources)} team resources in the view")
        
        return {
            'team_id': team_id,
            'molecule_id': molecule_id,
            'mixture_id': mixture_id,
            'invitation_id': invitation_id,
            'invitation_token': invitation_token,
            'resources': resources
        }
    except Exception as e:
        logger.error(f"Error creating test data: {e}")
        return {
            'error': str(e)
        }

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Implement proper team-based access control with RLS.")
    parser.add_argument("--test-only", action="store_true", help="Only test the existing team-based access control without applying the migration")
    parser.add_argument("--create-data", action="store_true", help="Create test data to validate team-based access control")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Create team tables if they don't exist
    logger.info("Checking teams structure...")
    team_structure = check_teams_structure()
    
    if not team_structure['teams_exists'] or not team_structure['team_members_exists']:
        logger.info("Creating team tables...")
        if not create_team_tables():
            logger.error("Failed to create team tables. Exiting.")
            sys.exit(1)
        
        # Refresh team structure
        team_structure = check_teams_structure()
    
    if not args.test_only:
        # Execute the migration
        logger.info("Executing migration to implement team-based access control...")
        success = execute_migration()
        
        if not success:
            logger.error("Migration failed. Exiting.")
            sys.exit(1)
    
    # Test team functions
    logger.info("Testing team functions...")
    function_tests = test_team_functions()
    
    # Test RLS policies
    logger.info("Testing RLS policies...")
    policy_tests = test_rls_policies()
    
    # Check team resources view
    logger.info("Checking team resources view...")
    view_check = check_team_view()
    
    # Create test data if requested
    test_data = None
    if args.create_data:
        logger.info("Creating test data...")
        test_data = create_test_data()
    
    # Save report
    report_file = save_report(
        team_structure,
        function_tests,
        policy_tests,
        view_check
    )
    
    if report_file:
        logger.info(f"Report saved to {report_file}")
    
    # Print summary
    structure_success = team_structure['teams_exists'] and team_structure['team_members_exists']
    function_success = all(v is True for k, v in function_tests.items() if k != 'error')
    policy_success = (
        len(policy_tests.get('teams_policies', [])) > 0 and
        len(policy_tests.get('team_members_policies', [])) > 0 and
        len(policy_tests.get('molecules_team_policies', [])) > 0 and
        len(policy_tests.get('mixtures_team_policies', [])) > 0
    )
    view_success = view_check.get('view_exists', False)
    
    overall_success = structure_success and function_success and policy_success and view_success
    
    if overall_success:
        logger.info("Team-based access control implementation successful!")
        print(f"\n✅ Team-based access control implementation successful!")
        print(f"   - Team tables structure: {'✅' if structure_success else '❌'}")
        print(f"   - Team functions: {'✅' if function_success else '❌'}")
        print(f"   - RLS policies: {'✅' if policy_success else '❌'}")
        print(f"   - Team resources view: {'✅' if view_success else '❌'}")
        
        if test_data and 'error' not in test_data:
            print(f"\n✅ Test data created successfully:")
            print(f"   - Team ID: {test_data['team_id']}")
            print(f"   - Test molecule ID: {test_data['molecule_id']}")
            print(f"   - Test mixture ID: {test_data['mixture_id']}")
            print(f"   - Test invitation token: {test_data['invitation_token']}")
        
        sys.exit(0)
    else:
        logger.warning("Team-based access control implementation has issues. Check the report for details.")
        print(f"\n⚠️ Team-based access control implementation has issues:")
        print(f"   - Team tables structure: {'✅' if structure_success else '❌'}")
        print(f"   - Team functions: {'✅' if function_success else '❌'}")
        print(f"   - RLS policies: {'✅' if policy_success else '❌'}")
        print(f"   - Team resources view: {'✅' if view_success else '❌'}")
        sys.exit(1)

if __name__ == "__main__":
    main()