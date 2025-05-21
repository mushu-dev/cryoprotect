#!/usr/bin/env python3
"""
CryoProtect - RDKit Service Tester

This script tests the RDKit service by sending requests to calculate properties
for various molecules and measuring the response time.
"""

import json
import time
import requests
import statistics
from datetime import datetime

# Configuration
RDKIT_URL = "http://localhost:5002"
TEST_MOLECULES = [
    {"name": "Ethanol", "smiles": "CCO"},
    {"name": "Glycerol", "smiles": "C(C(CO)O)O"},
    {"name": "DMSO", "smiles": "CS(=O)C"},
    {"name": "Sucrose", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)OC2C(C(C(C(O2)CO)O)O)O"}
]
NUM_ITERATIONS = 3

def test_health():
    """Test the health endpoint."""
    print(f"\n{'-' * 40}")
    print("Testing RDKit Service Health")
    print(f"{'-' * 40}")
    
    try:
        start_time = time.time()
        response = requests.get(f"{RDKIT_URL}/health", timeout=10)
        end_time = time.time()
        
        if response.status_code == 200:
            health_data = response.json()
            print(f"Status: {health_data.get('status', 'Unknown')}")
            print(f"RDKit Available: {health_data.get('rdkit_available', 'Unknown')}")
            print(f"RDKit Version: {health_data.get('rdkit_version', 'Unknown')}")
            print(f"RDKit Type: {health_data.get('rdkit_type', 'Unknown')}")
            print(f"Environment: {health_data.get('environment', 'Unknown')}")
            print(f"Response Time: {(end_time - start_time) * 1000:.2f} ms")
            return True
        else:
            print(f"Error: Service returned status code {response.status_code}")
            return False
    except Exception as e:
        print(f"Error: {str(e)}")
        return False

def test_property_calculation():
    """Test property calculation for various molecules."""
    print(f"\n{'-' * 40}")
    print("Testing Property Calculation")
    print(f"{'-' * 40}")
    
    results = {}
    
    for molecule in TEST_MOLECULES:
        name = molecule["name"]
        smiles = molecule["smiles"]
        
        print(f"\nTesting {name} ({smiles}):")
        
        response_times = []
        
        for i in range(NUM_ITERATIONS):
            try:
                start_time = time.time()
                response = requests.get(f"{RDKIT_URL}/calculate/{smiles}", timeout=10)
                end_time = time.time()
                response_time = (end_time - start_time) * 1000  # ms
                
                if response.status_code == 200:
                    data = response.json()
                    if i == 0:  # Only print properties on first iteration
                        print(f"  Properties:")
                        properties = data.get("properties", {})
                        for prop, value in properties.items():
                            print(f"    {prop}: {value}")
                    
                    print(f"  Iteration {i+1}: {response_time:.2f} ms")
                    response_times.append(response_time)
                else:
                    print(f"  Error: Service returned status code {response.status_code}")
            except Exception as e:
                print(f"  Error: {str(e)}")
                
        if response_times:
            avg_time = statistics.mean(response_times)
            min_time = min(response_times)
            max_time = max(response_times)
            
            results[name] = {
                "average": avg_time,
                "min": min_time,
                "max": max_time,
                "count": len(response_times)
            }
            
            print(f"  Summary: avg={avg_time:.2f} ms, min={min_time:.2f} ms, max={max_time:.2f} ms")
    
    return results

def run_performance_test():
    """Run a performance test with multiple concurrent requests."""
    print(f"\n{'-' * 40}")
    print("Running Performance Test")
    print(f"{'-' * 40}")
    
    import concurrent.futures
    
    # Create a list of 20 test requests (5 of each molecule)
    test_requests = []
    for _ in range(5):
        for molecule in TEST_MOLECULES:
            test_requests.append(molecule["smiles"])
    
    response_times = []
    errors = 0
    
    def fetch_properties(smiles):
        try:
            start_time = time.time()
            response = requests.get(f"{RDKIT_URL}/calculate/{smiles}", timeout=30)
            end_time = time.time()
            
            if response.status_code == 200:
                return (smiles, (end_time - start_time) * 1000, True)
            else:
                return (smiles, 0, False)
        except Exception:
            return (smiles, 0, False)
    
    print(f"Sending {len(test_requests)} concurrent requests...")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        start_time = time.time()
        futures = {executor.submit(fetch_properties, smiles): smiles for smiles in test_requests}
        
        for future in concurrent.futures.as_completed(futures):
            smiles, response_time, success = future.result()
            if success:
                response_times.append(response_time)
            else:
                errors += 1
        
        end_time = time.time()
    
    total_time = (end_time - start_time) * 1000
    
    if response_times:
        avg_time = statistics.mean(response_times)
        p95 = sorted(response_times)[int(0.95 * len(response_times))]
        
        print(f"\nPerformance Results:")
        print(f"  Total Time: {total_time:.2f} ms")
        print(f"  Average Response Time: {avg_time:.2f} ms")
        print(f"  95th Percentile: {p95:.2f} ms")
        print(f"  Requests: {len(test_requests)}")
        print(f"  Successful: {len(response_times)}")
        print(f"  Errors: {errors}")
        print(f"  Throughput: {len(response_times) / (total_time / 1000):.2f} req/sec")
        
        return {
            "total_time": total_time,
            "average_response_time": avg_time,
            "p95": p95,
            "requests": len(test_requests),
            "successful": len(response_times),
            "errors": errors,
            "throughput": len(response_times) / (total_time / 1000)
        }
    else:
        print("No successful responses")
        return None

def main():
    """Main function."""
    print(f"\n{'=' * 60}")
    print(f"CryoProtect RDKit Service Test - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'=' * 60}")
    
    # Test health
    if not test_health():
        print("\nRDKit service is not healthy. Aborting tests.")
        return
    
    # Test property calculation
    property_results = test_property_calculation()
    
    # Run performance test
    performance_results = run_performance_test()
    
    # Save results
    results = {
        "timestamp": datetime.now().isoformat(),
        "rdkit_url": RDKIT_URL,
        "property_results": property_results,
        "performance_results": performance_results
    }
    
    with open("rdkit_service_test_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'=' * 60}")
    print("Tests completed. Results saved to rdkit_service_test_results.json")
    print(f"{'=' * 60}")

if __name__ == "__main__":
    main()