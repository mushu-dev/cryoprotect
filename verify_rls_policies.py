import requests
import json
import time
import statistics
from typing import Dict, List, Any, Optional, Tuple
import uuid

# Configuration
SUPABASE_URL = "https://tsdlmynydfuypiugmkev.supabase.co"
ANON_KEY = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InRzZGxteW55ZGZ1eXBpdWdta2V2Iiwicm9sZSI6ImFub24iLCJpYXQiOjE3NDQ3ODAxODMsImV4cCI6MjA2MDM1NjE4M30.T6udmD-3lhlTy4bY2p0y2lX5-11yvyn425PWrlnIPLU"

# Test users - these should be created in the Supabase Auth system
# For testing purposes, we'll simulate different roles
TEST_USERS = {
    "anon": {
        "role": "anon",
        "token": ANON_KEY
    },
    # We'll simulate authenticated users with different roles
    "authenticated": {
        "role": "authenticated",
        "token": None  # Will be populated after login
    },
    "service_role": {
        "role": "service_role",
        "token": None  # This would be the service_role key in a real scenario
    }
}

# Tables to test
TABLES_TO_TEST = [
    "molecules", 
    "mixtures", 
    "mixture_components", 
    "experiments", 
    "experiment_properties",
    "molecular_properties", 
    "predictions", 
    "property_types",
    "calculation_methods",
    "proteins",
    "molecule_proteins",
    "molecule_experiments",
    "scientific_data_audit"
]

class SupabaseClient:
    def __init__(self, url: str, key: str):
        self.url = url
        self.key = key
        self.headers = {
            "apikey": key,
            "Authorization": f"Bearer {key}",
            "Content-Type": "application/json"
        }
    
    def select(self, table: str, limit: int = 10) -> Dict:
        """Select data from a table"""
        response = requests.get(
            f"{self.url}/rest/v1/{table}?select=*&limit={limit}",
            headers=self.headers
        )
        return response.json() if response.status_code == 200 else {"error": response.text, "status": response.status_code}
    
    def insert(self, table: str, data: Dict) -> Dict:
        """Insert data into a table"""
        response = requests.post(
            f"{self.url}/rest/v1/{table}",
            headers=self.headers,
            json=data
        )
        return response.json() if response.status_code == 201 else {"error": response.text, "status": response.status_code}
    
    def update(self, table: str, id_field: str, id_value: str, data: Dict) -> Dict:
        """Update data in a table"""
        response = requests.patch(
            f"{self.url}/rest/v1/{table}?{id_field}=eq.{id_value}",
            headers=self.headers,
            json=data
        )
        return {"success": True} if response.status_code == 204 else {"error": response.text, "status": response.status_code}
    
    def delete(self, table: str, id_field: str, id_value: str) -> Dict:
        """Delete data from a table"""
        response = requests.delete(
            f"{self.url}/rest/v1/{table}?{id_field}=eq.{id_value}",
            headers=self.headers
        )
        return {"success": True} if response.status_code == 204 else {"error": response.text, "status": response.status_code}
    
    def execute_sql(self, query: str) -> Dict:
        """Execute a raw SQL query using the REST API"""
        # This is a simplified version - in a real scenario, you'd use a server-side function
        # or the Supabase client's rpc method
        response = requests.post(
            f"{self.url}/rest/v1/rpc/exec_sql",
            headers=self.headers,
            json={"query": query}
        )
        return response.json() if response.status_code == 200 else {"error": response.text, "status": response.status_code}
    
    def measure_query_performance(self, table: str, iterations: int = 5) -> Dict:
        """Measure query performance"""
        times = []
        for _ in range(iterations):
            start_time = time.time()
            self.select(table, limit=100)
            end_time = time.time()
            times.append(end_time - start_time)
        
        return {
            "min": min(times),
            "max": max(times),
            "avg": sum(times) / len(times),
            "median": statistics.median(times),
            "std_dev": statistics.stdev(times) if len(times) > 1 else 0
        }

class RLSVerifier:
    def __init__(self):
        self.results = {
            "rls_effectiveness": {},
            "access_patterns": {},
            "performance": {},
            "data_relationships": {}
        }
        self.clients = {
            role: SupabaseClient(SUPABASE_URL, user_info["token"])
            for role, user_info in TEST_USERS.items()
            if user_info["token"] is not None
        }
    
    def test_rls_effectiveness(self) -> Dict:
        """Test RLS effectiveness on all tables"""
        results = {}
        
        for table in TABLES_TO_TEST:
            table_results = {}
            
            # Test with anonymous user
            if "anon" in self.clients:
                anon_client = self.clients["anon"]
                
                # Test SELECT
                select_result = anon_client.select(table)
                table_results["anon_select"] = {
                    "success": "error" not in select_result,
                    "data_count": len(select_result) if "error" not in select_result else 0,
                    "error": select_result.get("error", None)
                }
                
                # Test INSERT (should fail for most tables with RLS)
                insert_data = {"name": f"Test {uuid.uuid4()}", "is_public": True}
                insert_result = anon_client.insert(table, insert_data)
                table_results["anon_insert"] = {
                    "success": "error" not in insert_result,
                    "error": insert_result.get("error", None)
                }
            
            # Test with service role (if available)
            if "service_role" in self.clients:
                service_client = self.clients["service_role"]
                
                # Test SELECT
                select_result = service_client.select(table)
                table_results["service_role_select"] = {
                    "success": "error" not in select_result,
                    "data_count": len(select_result) if "error" not in select_result else 0,
                    "error": select_result.get("error", None)
                }
            
            results[table] = table_results
        
        self.results["rls_effectiveness"] = results
        return results
    
    def test_access_patterns(self) -> Dict:
        """Test access patterns with different user roles"""
        results = {}
        
        # Test public vs private data access
        for table in ["molecules", "mixtures", "proteins"]:
            table_results = {}
            
            if "anon" in self.clients:
                anon_client = self.clients["anon"]
                
                # Get public records
                public_records = anon_client.select(table)
                if "error" not in public_records:
                    public_count = sum(1 for record in public_records if record.get("is_public", False))
                    table_results["public_records_accessible"] = public_count
                
                # Try to access a specific non-public record (if we know one)
                # This would require knowing a specific ID of a non-public record
                # For now, we'll just check if any non-public records are returned
                non_public_count = sum(1 for record in public_records if not record.get("is_public", True))
                table_results["non_public_records_accessible"] = non_public_count
            
            results[table] = table_results
        
        # Test relationship-based access
        relationship_tests = [
            ("molecules", "molecular_properties", "molecule_id"),
            ("mixtures", "mixture_components", "mixture_id"),
            ("experiments", "experiment_properties", "experiment_id")
        ]
        
        for parent_table, child_table, foreign_key in relationship_tests:
            test_key = f"{parent_table}_{child_table}"
            table_results = {}
            
            if "anon" in self.clients:
                anon_client = self.clients["anon"]
                
                # Get accessible parent records
                parent_records = anon_client.select(parent_table)
                if "error" not in parent_records and len(parent_records) > 0:
                    parent_id = parent_records[0]["id"]
                    
                    # Try to access related child records
                    child_records = requests.get(
                        f"{SUPABASE_URL}/rest/v1/{child_table}?{foreign_key}=eq.{parent_id}&select=*",
                        headers=anon_client.headers
                    ).json()
                    
                    table_results["related_records_accessible"] = {
                        "parent_id": parent_id,
                        "child_records_count": len(child_records) if isinstance(child_records, list) else 0
                    }
            
            results[test_key] = table_results
        
        self.results["access_patterns"] = results
        return results
    
    def test_query_performance(self) -> Dict:
        """Measure query performance with RLS enabled"""
        results = {}
        
        for table in TABLES_TO_TEST:
            if "anon" in self.clients:
                results[table] = self.clients["anon"].measure_query_performance(table)
        
        self.results["performance"] = results
        return results
    
    def validate_data_relationships(self) -> Dict:
        """Validate scientific data relationships (referential integrity, foreign keys)"""
        results = {}
        
        # Define relationships to test
        relationships = [
            ("mixture_components", "mixture_id", "mixtures", "id"),
            ("mixture_components", "molecule_id", "molecules", "id"),
            ("molecular_properties", "molecule_id", "molecules", "id"),
            ("molecular_properties", "property_type_id", "property_types", "id"),
            ("experiments", "molecule_id", "molecules", "id"),
            ("experiments", "mixture_id", "mixtures", "id"),
            ("experiments", "property_type_id", "property_types", "id"),
            ("experiment_properties", "experiment_id", "experiments", "id"),
            ("experiment_properties", "property_type_id", "property_types", "id"),
            ("predictions", "molecule_id", "molecules", "id"),
            ("predictions", "mixture_id", "mixtures", "id"),
            ("predictions", "property_type_id", "property_types", "id"),
            ("predictions", "calculation_method_id", "calculation_methods", "id"),
            ("molecule_proteins", "molecule_id", "molecules", "id"),
            ("molecule_proteins", "protein_id", "proteins", "id"),
            ("molecule_experiments", "molecule_id", "molecules", "id"),
            ("molecule_experiments", "experiment_id", "experiments", "id")
        ]
        
        if "service_role" in self.clients:
            service_client = self.clients["service_role"]
            
            for child_table, fk_column, parent_table, pk_column in relationships:
                # This would ideally be a SQL query to check referential integrity
                # For now, we'll just check if we can access the related records
                child_records = service_client.select(child_table, limit=5)
                
                if "error" not in child_records and len(child_records) > 0:
                    valid_relationships = 0
                    invalid_relationships = 0
                    
                    for child in child_records:
                        if fk_column in child and child[fk_column] is not None:
                            fk_value = child[fk_column]
                            parent_record = requests.get(
                                f"{SUPABASE_URL}/rest/v1/{parent_table}?{pk_column}=eq.{fk_value}&select=*",
                                headers=service_client.headers
                            ).json()
                            
                            if isinstance(parent_record, list) and len(parent_record) > 0:
                                valid_relationships += 1
                            else:
                                invalid_relationships += 1
                    
                    results[f"{child_table}_{fk_column}"] = {
                        "valid_relationships": valid_relationships,
                        "invalid_relationships": invalid_relationships,
                        "integrity_percentage": (valid_relationships / (valid_relationships + invalid_relationships) * 100) if (valid_relationships + invalid_relationships) > 0 else 0
                    }
        
        self.results["data_relationships"] = results
        return results
    
    def run_all_tests(self) -> Dict:
        """Run all verification tests"""
        print("Testing RLS effectiveness...")
        self.test_rls_effectiveness()
        
        print("Testing access patterns...")
        self.test_access_patterns()
        
        print("Testing query performance...")
        self.test_query_performance()
        
        print("Validating data relationships...")
        self.validate_data_relationships()
        
        return self.results
    
    def generate_report(self, output_format: str = "json") -> str:
        """Generate a comprehensive verification report"""
        if output_format == "json":
            return json.dumps(self.results, indent=2)
        elif output_format == "markdown":
            md = "# CryoProtect Database Verification Report\n\n"
            
            # RLS Effectiveness
            md += "## 1. RLS Effectiveness\n\n"
            for table, results in self.results["rls_effectiveness"].items():
                md += f"### {table}\n\n"
                md += "| Test | Success | Data Count | Error |\n"
                md += "|------|---------|------------|-------|\n"
                
                for test, test_results in results.items():
                    success = test_results.get("success", False)
                    data_count = test_results.get("data_count", "N/A")
                    error = test_results.get("error", "None")
                    md += f"| {test} | {success} | {data_count} | {error} |\n"
                
                md += "\n"
            
            # Access Patterns
            md += "## 2. Access Patterns\n\n"
            for test, results in self.results["access_patterns"].items():
                md += f"### {test}\n\n"
                md += "```json\n"
                md += json.dumps(results, indent=2)
                md += "\n```\n\n"
            
            # Performance
            md += "## 3. Query Performance (seconds)\n\n"
            md += "| Table | Min | Max | Average | Median | Std Dev |\n"
            md += "|-------|-----|-----|---------|--------|--------|\n"
            
            for table, metrics in self.results["performance"].items():
                min_time = metrics.get("min", "N/A")
                max_time = metrics.get("max", "N/A")
                avg_time = metrics.get("avg", "N/A")
                median_time = metrics.get("median", "N/A")
                std_dev = metrics.get("std_dev", "N/A")
                
                md += f"| {table} | {min_time:.4f} | {max_time:.4f} | {avg_time:.4f} | {median_time:.4f} | {std_dev:.4f} |\n"
            
            md += "\n"
            
            # Data Relationships
            md += "## 4. Data Relationships\n\n"
            md += "| Relationship | Valid | Invalid | Integrity % |\n"
            md += "|--------------|-------|---------|------------|\n"
            
            for relationship, metrics in self.results["data_relationships"].items():
                valid = metrics.get("valid_relationships", 0)
                invalid = metrics.get("invalid_relationships", 0)
                integrity = metrics.get("integrity_percentage", 0)
                
                md += f"| {relationship} | {valid} | {invalid} | {integrity:.2f}% |\n"
            
            return md
        else:
            return "Unsupported output format"

if __name__ == "__main__":
    verifier = RLSVerifier()
    results = verifier.run_all_tests()
    
    # Generate reports
    json_report = verifier.generate_report("json")
    md_report = verifier.generate_report("markdown")
    
    # Save reports
    with open("rls_verification_report.json", "w") as f:
        f.write(json_report)
    
    with open("RLS_Verification_Report.md", "w") as f:
        f.write(md_report)
    
    print("Verification complete. Reports saved to rls_verification_report.json and RLS_Verification_Report.md")