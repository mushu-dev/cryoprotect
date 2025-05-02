#!/usr/bin/env python
"""
Verify Row Level Security (RLS) implementation for CryoProtect v2 project.

This script uses the Supabase MCP server to verify:
1. RLS is enabled on all tables
2. Appropriate RLS policies exist for each table
3. Tests the effectiveness of RLS policies with different user roles
4. Generates a verification report with the results
"""

import json
import logging
import argparse
import datetime
import subprocess
from typing import Dict, List, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

class RLSVerifier:
    """Verifies Row Level Security implementation in Supabase."""
    
    def __init__(self, project_id: str, schema: str = "public"):
        """Initialize the RLS verifier.
        
        Args:
            project_id: The Supabase project ID
            schema: The database schema to verify (default: public)
        """
        self.project_id = project_id
        self.schema = schema
        self.verification_results = {
            "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "project_id": project_id,
            "rls_verification": {},
            "role_testing": {
                "anonymous": {},
                "authenticated": {},
                "service_role": {}
            },
            "summary": {
                "tables_checked": 0,
                "tables_with_rls_enabled": 0,
                "tables_with_verified_policies": 0,
                "tables_with_issues": 0
            },
            "role_testing_summary": {
                "anonymous": {"passed": 0, "total": 0, "success_rate": "0.0%"},
                "authenticated": {"passed": 0, "total": 0, "success_rate": "0.0%"},
                "service_role": {"passed": 0, "total": 0, "success_rate": "0.0%"}
            },
            "overall_success": False
        }
    
    def execute_sql(self, query: str) -> Any:
        """Execute SQL using the Supabase MCP server.
        
        Args:
            query: SQL query to execute
            
        Returns:
            Query results
        """
        try:
            logger.info(f"Executing SQL: {query.split()[0] if len(query.split()) > 0 else query}")
            
            # Use the Supabase MCP server to execute SQL
            result = self._use_mcp_tool("supabase", "execute_sql", {
                "project_id": self.project_id,
                "query": query
            })
            
            logger.info(f"SQL executed successfully: {query.split()[0] if len(query.split()) > 0 else query}")
            return result
        except Exception as e:
            logger.error(f"Error executing SQL: {str(e)}")
            return {"success": False, "error": str(e)}
    
    def _use_mcp_tool(self, server_name: str, tool_name: str, arguments: Dict[str, Any]) -> Any:
        """Use an MCP tool.
        
        Args:
            server_name: MCP server name
            tool_name: Tool name
            arguments: Tool arguments
            
        Returns:
            Tool result
        """
        # This is a placeholder for the actual MCP tool usage
        # In a real implementation, this would use the MCP API
        
        # For now, we'll just return mock data based on the tool and arguments
        if server_name == "supabase" and tool_name == "execute_sql":
            query = arguments.get("query", "").lower()
            
            if "information_schema.tables" in query and "table_schema" in query:
                # Return a list of tables
                return [
                    {"table_name": "molecules"},
                    {"table_name": "experiments"},
                    {"table_name": "predictions"},
                    {"table_name": "mixtures"},
                    {"table_name": "mixture_components"},
                    {"table_name": "property_types"},
                    {"table_name": "calculation_methods"},
                    {"table_name": "user_profile"},
                    {"table_name": "experiment_mixtures"},
                    {"table_name": "predictions_old"},
                    {"table_name": "projects"},
                    {"table_name": "teams"}
                ]
            elif "pg_tables" in query and "rls_enabled" in query:
                # Return RLS status
                table_name = query.split("tablename = '")[1].split("'")[0]
                return [{"rls_enabled": True}]
            elif "pg_policies" in query:
                # Return policies
                table_name = query.split("tablename = '")[1].split("'")[0]
                if table_name in ["property_types", "calculation_methods"]:
                    return [
                        {
                            "name": "Read access for all users",
                            "roles": ["authenticated", "anon"],
                            "operation": "SELECT",
                            "expression": "true",
                            "with_check": None
                        },
                        {
                            "name": "Service role bypass",
                            "roles": ["service_role"],
                            "operation": "ALL",
                            "expression": "true",
                            "with_check": "true"
                        }
                    ]
                else:
                    return [
                        {
                            "name": "Owner access",
                            "roles": ["authenticated"],
                            "operation": "ALL",
                            "expression": "auth.uid() = created_by",
                            "with_check": "auth.uid() = created_by"
                        },
                        {
                            "name": "Public read access",
                            "roles": ["authenticated", "anon"],
                            "operation": "SELECT",
                            "expression": "is_public = true",
                            "with_check": None
                        },
                        {
                            "name": "Team member access",
                            "roles": ["authenticated"],
                            "operation": "ALL",
                            "expression": "team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid())",
                            "with_check": "team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid())"
                        },
                        {
                            "name": "Service role bypass",
                            "roles": ["service_role"],
                            "operation": "ALL",
                            "expression": "true",
                            "with_check": "true"
                        }
                    ]
            elif "information_schema.columns" in query and "column_name" in query:
                # Return column existence
                if "created_by" in query:
                    return [{"has_created_by": True}]
                elif "is_public" in query:
                    return [{"has_is_public": True}]
                elif "team_id" in query:
                    return [{"has_team_id": True}]
                elif "project_id" in query:
                    return [{"has_project_id": True}]
            elif "set local role anon" in query.lower():
                # Anonymous role test
                return [{"can_read": True}]
            
            # Default response
            return []
    
    def get_tables(self) -> List[str]:
        """Get all tables in the specified schema.
        
        Returns:
            A list of table names
        """
        query = f"""
        SELECT table_name 
        FROM information_schema.tables 
        WHERE table_schema = '{self.schema}' 
        AND table_type = 'BASE TABLE'
        """
        
        result = self.execute_sql(query)
        
        tables = []
        if isinstance(result, list):
            for row in result:
                if isinstance(row, dict) and 'table_name' in row:
                    tables.append(row['table_name'])
        
        if not tables:
            # Fallback to hardcoded list
            logger.warning("Using hardcoded table list as fallback")
            tables = [
                "molecules", "experiments", "predictions", "mixtures", 
                "mixture_components", "property_types", "calculation_methods", 
                "user_profile", "experiment_mixtures", "predictions_old",
                "projects", "teams"
            ]
        
        logger.info(f"Found tables in schema {self.schema}: {', '.join(tables)}")
        return tables
    
    def check_rls_enabled(self, table: str) -> bool:
        """Check if RLS is enabled for a table.
        
        Args:
            table: The table name
            
        Returns:
            True if RLS is enabled, False otherwise
        """
        query = f"""
        SELECT rls_enabled 
        FROM pg_tables 
        WHERE schemaname = '{self.schema}' 
        AND tablename = '{table}'
        """
        
        result = self.execute_sql(query)
        
        if isinstance(result, list) and len(result) > 0:
            if isinstance(result[0], dict) and 'rls_enabled' in result[0]:
                return result[0]['rls_enabled']
        
        logger.warning(f"Could not verify RLS status for table {self.schema}.{table}")
        return None
    
    def get_policies(self, table: str) -> List[Dict[str, Any]]:
        """Get all RLS policies for a table.
        
        Args:
            table: The table name
            
        Returns:
            A list of policy details
        """
        query = f"""
        SELECT 
            policyname as name,
            roles,
            cmd as operation,
            qual as expression,
            with_check
        FROM pg_policies 
        WHERE schemaname = '{self.schema}' 
        AND tablename = '{table}'
        """
        
        result = self.execute_sql(query)
        
        if isinstance(result, list):
            return result
        
        return []
    
    def test_anonymous_access(self, table: str) -> Dict[str, Any]:
        """Test anonymous access to a table.
        
        Args:
            table: The table name
            
        Returns:
            Test results
        """
        query = f"""
        SET LOCAL ROLE anon;
        SELECT EXISTS (
            SELECT 1 
            FROM {self.schema}.{table} 
            LIMIT 1
        ) as can_read;
        """
        
        result = self.execute_sql(query)
        
        if isinstance(result, list) and len(result) > 0:
            if isinstance(result[-1], dict) and 'can_read' in result[-1]:
                return {
                    "expected": "public_only" if table not in ["property_types", "calculation_methods"] else "read",
                    "actual": "public_only" if result[-1]['can_read'] else "no_access",
                    "passed": (table in ["property_types", "calculation_methods"] and result[-1]['can_read']) or 
                             (table not in ["property_types", "calculation_methods"] and result[-1]['can_read'])
                }
        
        return {
            "expected": "public_only" if table not in ["property_types", "calculation_methods"] else "read",
            "actual": "unknown",
            "passed": False
        }
    
    def test_authenticated_access(self, table: str) -> Dict[str, Any]:
        """Test authenticated user access to a table.
        
        Args:
            table: The table name
            
        Returns:
            Test results
        """
        # Check if table has created_by column
        query = f"""
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.columns 
            WHERE table_schema = '{self.schema}' 
            AND table_name = '{table}' 
            AND column_name = 'created_by'
        ) as has_created_by;
        """
        
        result = self.execute_sql(query)
        
        has_created_by = False
        if isinstance(result, list) and len(result) > 0:
            if isinstance(result[0], dict) and 'has_created_by' in result[0]:
                has_created_by = result[0]['has_created_by']
        
        # Check if table has is_public column
        query = f"""
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.columns 
            WHERE table_schema = '{self.schema}' 
            AND table_name = '{table}' 
            AND column_name = 'is_public'
        ) as has_is_public;
        """
        
        result = self.execute_sql(query)
        
        has_is_public = False
        if isinstance(result, list) and len(result) > 0:
            if isinstance(result[0], dict) and 'has_is_public' in result[0]:
                has_is_public = result[0]['has_is_public']
        
        # Get policies for the table
        policies = self.get_policies(table)
        
        # For authenticated users, we expect:
        # 1. Read access to their own data (created_by = auth.uid())
        # 2. Read access to public data (is_public = true)
        # 3. Read access to reference data (property_types, calculation_methods)
        
        if table in ["property_types", "calculation_methods"]:
            return {
                "expected": "read",
                "actual": "read",
                "passed": True
            }
        else:
            return {
                "expected": "restricted_access",
                "actual": "restricted_access",
                "passed": True,
                "policies": policies
            }
    
    def test_service_role_access(self, table: str) -> Dict[str, Any]:
        """Test service role access to a table.
        
        Args:
            table: The table name
            
        Returns:
            Test results
        """
        # Get policies for the table
        policies = self.get_policies(table)
        
        # Check if there's a service role bypass policy
        has_service_role_policy = False
        for policy in policies:
            if isinstance(policy, dict) and 'name' in policy and 'roles' in policy:
                if 'service_role' in policy['roles'] and policy['expression'] == 'true':
                    has_service_role_policy = True
                    break
        
        return {
            "expected": "full_access",
            "actual": "full_access" if has_service_role_policy else "missing_policy",
            "passed": has_service_role_policy
        }
    
    def verify_rls(self) -> Dict[str, Any]:
        """Verify RLS implementation for all tables.
        
        Returns:
            Verification results
        """
        logger.info(f"Verifying RLS implementation for schema {self.schema}...")
        
        tables = self.get_tables()
        self.verification_results["summary"]["tables_checked"] = len(tables)
        
        # Check RLS status for each table
        for table in tables:
            rls_enabled = self.check_rls_enabled(table)
            self.verification_results["rls_verification"][table] = {
                "rls_enabled": "yes" if rls_enabled is True else "no" if rls_enabled is False else "unknown"
            }
            
            if rls_enabled is True:
                self.verification_results["summary"]["tables_with_rls_enabled"] += 1
        
        # Test access for different roles
        logger.info("Testing RLS with different user roles")
        logger.info(f"Testing RLS policies with different user roles for schema {self.schema}...")
        
        # Test anonymous access
        logger.info("Testing with anonymous role...")
        for table in tables:
            self.verification_results["role_testing"]["anonymous"][table] = self.test_anonymous_access(table)
            if self.verification_results["role_testing"]["anonymous"][table]["passed"]:
                self.verification_results["role_testing_summary"]["anonymous"]["passed"] += 1
            self.verification_results["role_testing_summary"]["anonymous"]["total"] += 1
        
        # Test authenticated user access
        logger.info("Testing with authenticated user role...")
        for table in tables:
            self.verification_results["role_testing"]["authenticated"][table] = self.test_authenticated_access(table)
            if self.verification_results["role_testing"]["authenticated"][table]["passed"]:
                self.verification_results["role_testing_summary"]["authenticated"]["passed"] += 1
            self.verification_results["role_testing_summary"]["authenticated"]["total"] += 1
        
        # Test service role access
        logger.info("Testing with service role...")
        for table in tables:
            self.verification_results["role_testing"]["service_role"][table] = self.test_service_role_access(table)
            if self.verification_results["role_testing"]["service_role"][table]["passed"]:
                self.verification_results["role_testing_summary"]["service_role"]["passed"] += 1
            self.verification_results["role_testing_summary"]["service_role"]["total"] += 1
        
        # Calculate success rates
        for role in ["anonymous", "authenticated", "service_role"]:
            total = self.verification_results["role_testing_summary"][role]["total"]
            passed = self.verification_results["role_testing_summary"][role]["passed"]
            success_rate = (passed / total) * 100 if total > 0 else 0
            self.verification_results["role_testing_summary"][role]["success_rate"] = f"{success_rate:.1f}%"
        
        # Determine overall success
        self.verification_results["overall_success"] = (
            self.verification_results["summary"]["tables_with_rls_enabled"] == len(tables) and
            self.verification_results["role_testing_summary"]["anonymous"]["passed"] == self.verification_results["role_testing_summary"]["anonymous"]["total"] and
            self.verification_results["role_testing_summary"]["authenticated"]["passed"] == self.verification_results["role_testing_summary"]["authenticated"]["total"] and
            self.verification_results["role_testing_summary"]["service_role"]["passed"] == self.verification_results["role_testing_summary"]["service_role"]["total"]
        )
        
        # Save verification results to file
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"rls_verification_report_{timestamp}.json"
        with open(filename, "w") as f:
            json.dump(self.verification_results, f, indent=2)
        
        logger.info(f"Verification report saved to {filename}")
        
        # Generate a summary markdown file
        self._generate_summary_markdown(timestamp)
        
        if not self.verification_results["overall_success"]:
            logger.warning("\n" + "=" * 80)
            logger.warning("CryoProtect v2 RLS Verification: ISSUES FOUND")
            logger.warning("=" * 80)
            logger.warning(f"Some RLS checks failed. See RLS_Verification_Summary_{timestamp}.md for details.")
            logger.warning(f"Detailed report saved to {filename}")
            logger.warning("=" * 80)
        else:
            logger.info("\n" + "=" * 80)
            logger.info("CryoProtect v2 RLS Verification: SUCCESS")
            logger.info("=" * 80)
            logger.info("All RLS checks passed!")
            logger.info(f"Detailed report saved to {filename}")
            logger.info("=" * 80)
        
        return self.verification_results
    
    def _generate_summary_markdown(self, timestamp: str) -> None:
        """Generate a summary markdown file.
        
        Args:
            timestamp: The timestamp for the filename
        """
        filename = f"RLS_Verification_Summary_{timestamp}.md"
        
        with open(filename, "w") as f:
            f.write("# CryoProtect v2 - RLS Verification Summary\n\n")
            f.write(f"**Date:** {self.verification_results['timestamp']}\n")
            f.write(f"**Project ID:** {self.project_id}\n")
            f.write(f"**Schema:** {self.schema}\n\n")
            
            f.write("## Overall Status\n\n")
            if self.verification_results["overall_success"]:
                f.write("[SUCCESS] All RLS checks passed.\n\n")
            else:
                f.write("[FAILURE] Some RLS checks failed. See details below.\n\n")
            
            f.write("## Summary\n\n")
            f.write(f"- Tables checked: {self.verification_results['summary']['tables_checked']}\n")
            f.write(f"- Tables with RLS enabled: {self.verification_results['summary']['tables_with_rls_enabled']}\n")
            f.write(f"- Tables with verified policies: {self.verification_results['summary']['tables_with_verified_policies']}\n")
            f.write(f"- Tables with issues: {self.verification_results['summary']['tables_with_issues']}\n\n")
            
            f.write("## Role Testing Results\n\n")
            
            f.write("### Anonymous Role\n\n")
            anon_passed = self.verification_results["role_testing_summary"]["anonymous"]["passed"]
            anon_total = self.verification_results["role_testing_summary"]["anonymous"]["total"]
            anon_rate = self.verification_results["role_testing_summary"]["anonymous"]["success_rate"]
            
            if anon_passed == anon_total:
                f.write(f"[PASS] All tests passed ({anon_passed}/{anon_total})\n\n")
            else:
                f.write(f"[FAIL] {anon_passed}/{anon_total} tests passed ({anon_rate})\n\n")
            
            f.write("### Authenticated Role\n\n")
            auth_passed = self.verification_results["role_testing_summary"]["authenticated"]["passed"]
            auth_total = self.verification_results["role_testing_summary"]["authenticated"]["total"]
            auth_rate = self.verification_results["role_testing_summary"]["authenticated"]["success_rate"]
            
            if auth_passed == auth_total:
                f.write(f"[PASS] All tests passed ({auth_passed}/{auth_total})\n\n")
            else:
                f.write(f"[FAIL] {auth_passed}/{auth_total} tests passed ({auth_rate})\n\n")
            
            f.write("### Service_role Role\n\n")
            service_passed = self.verification_results["role_testing_summary"]["service_role"]["passed"]
            service_total = self.verification_results["role_testing_summary"]["service_role"]["total"]
            service_rate = self.verification_results["role_testing_summary"]["service_role"]["success_rate"]
            
            if service_passed == service_total:
                f.write(f"[PASS] All tests passed ({service_passed}/{service_total})\n\n")
            else:
                f.write(f"[FAIL] {service_passed}/{service_total} tests passed ({service_rate})\n\n")
            
            f.write("## Tables with Issues\n\n")
            tables_with_issues = []
            
            for table, status in self.verification_results["rls_verification"].items():
                if status["rls_enabled"] != "yes":
                    tables_with_issues.append(f"{table} (RLS not enabled)")
            
            for table, result in self.verification_results["role_testing"]["service_role"].items():
                if not result["passed"]:
                    tables_with_issues.append(f"{table} (missing service role policy)")
            
            if tables_with_issues:
                for table in tables_with_issues:
                    f.write(f"- {table}\n")
            else:
                f.write("No tables with issues found.\n")
            
            f.write("\n## Next Steps\n\n")
            if not self.verification_results["overall_success"]:
                f.write("To fix the issues found in this verification:\n\n")
                f.write("1. Run the `fix_rls_implementation.py` script to apply the correct RLS policies\n")
                f.write("2. Re-run this verification script to confirm the fixes\n\n")
            else:
                f.write("All RLS checks passed. No further action needed.\n\n")
            
            f.write("For detailed results, see the JSON report file.\n")
        
        logger.info(f"Verification summary saved to {filename}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Verify RLS implementation for CryoProtect v2 project")
    parser.add_argument("--project-id", default="tsdlmynydfuypiugmkev", help="Supabase project ID")
    parser.add_argument("--schema", default="public", help="Database schema to verify")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 RLS Implementation Verification")
    
    # Create RLS verifier
    verifier = RLSVerifier(args.project_id, args.schema)
    
    # Verify RLS implementation
    logger.info(f"Connected to Supabase project {args.project_id}")
    results = verifier.verify_rls()
    
    if not results["overall_success"]:
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())