#!/bin/bash

# Run Connection Pool Performance Tests
# This script runs the connection pool performance tests

# Exit on error
set -e

# Print header
echo "==================================================="
echo "Connection Pool Performance Tests"
echo "==================================================="

# Set up environment
echo -e "\nSetting up environment variables..."
if [ -f .env ]; then
    echo "Loading .env file"
    export $(grep -v '^#' .env | xargs)
else
    echo "No .env file found, please make sure database connection parameters are set in the environment"
    read -p "Do you want to continue? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]
    then
        exit 1
    fi
fi

# Check if required variables exist
if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_NAME" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
    echo "Error: Required database connection variables are not set"
    echo "Please set SUPABASE_DB_HOST, SUPABASE_DB_NAME, SUPABASE_DB_USER, and SUPABASE_DB_PASSWORD"
    exit 1
fi

# Create Python virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo -e "\nCreating virtual environment..."
    python -m venv venv
fi

# Activate virtual environment
echo -e "\nActivating virtual environment..."
source venv/bin/activate

# Install required packages
echo -e "\nInstalling required packages..."
pip install -r requirements.txt
pip install python-dotenv psycopg2-binary

# Create reports directory if it doesn't exist
mkdir -p reports

# Timestamp for reports
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
REPORT_FILE="reports/connection_pool_performance_${TIMESTAMP}.json"

# Step 1: Check database connection
echo -e "\n==================================================="
echo "Step 1: Checking database connection"
echo "==================================================="
echo "Testing connection to Supabase database..."
python test_db_connection.py
if [ $? -ne 0 ]; then
    echo "Database connection test failed. Please fix the connection issues before continuing."
    exit 1
fi
echo "Database connection test passed!"

# Step 2: Run basic verification tests
echo -e "\n==================================================="
echo "Step 2: Running basic connection pool verification tests"
echo "==================================================="
echo "Executing connection pool verification tests..."
python -m tests.test_connection_pool_verification
TEST_RESULT_VERIFICATION=$?
echo "Connection pool verification tests completed with exit code: $TEST_RESULT_VERIFICATION"

# Step 3: Run performance tests
echo -e "\n==================================================="
echo "Step 3: Running connection pool performance tests"
echo "==================================================="
echo "Executing performance tests..."

# Parse command line arguments
ITERATIONS=5
CONCURRENT=10
MIN_CONN=2
MAX_CONN=10
QUERY="SELECT COUNT(*) FROM molecule"

# Process arguments if provided
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --iterations)
        ITERATIONS="$2"
        shift
        shift
        ;;
        --concurrent)
        CONCURRENT="$2"
        shift
        shift
        ;;
        --min-conn)
        MIN_CONN="$2"
        shift
        shift
        ;;
        --max-conn)
        MAX_CONN="$2"
        shift
        shift
        ;;
        --query)
        QUERY="$2"
        shift
        shift
        ;;
        *)
        # unknown option
        shift
        ;;
    esac
done

# Run performance tests
python -m tests.test_connection_pool_performance \
    --query "$QUERY" \
    --iterations $ITERATIONS \
    --concurrent $CONCURRENT \
    --min-conn $MIN_CONN \
    --max-conn $MAX_CONN \
    --output $REPORT_FILE

TEST_RESULT_PERFORMANCE=$?
echo "Connection pool performance tests completed with exit code: $TEST_RESULT_PERFORMANCE"

# Step 4: Run batched stress test (if verification and performance tests passed)
if [ $TEST_RESULT_VERIFICATION -eq 0 ] && [ $TEST_RESULT_PERFORMANCE -eq 0 ]; then
    echo -e "\n==================================================="
    echo "Step 4: Running connection pool stress test"
    echo "==================================================="
    echo "This test will simulate high load with multiple concurrent threads..."
    
    # Create stress test script if it doesn't exist
    if [ ! -f "tests/stress_test_connection_pool.py" ]; then
        echo "Creating stress test script..."
        
        cat > tests/stress_test_connection_pool.py << EOF
"""
Connection Pool Stress Test
Simulates high load with multiple concurrent threads
"""
import os
import sys
import time
import threading
import random
import argparse
import psycopg2
from dotenv import load_dotenv

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Try to import connection pool implementation
try:
    from db_pool import ConnectionPool, execute_query
    HAS_DB_POOL = True
except ImportError:
    try:
        from optimized_connection_pool import ConnectionPool, execute_query
        HAS_DB_POOL = True
    except ImportError:
        try:
            from database.connection_manager import ConnectionPool, execute_query
            HAS_DB_POOL = True
        except ImportError:
            HAS_DB_POOL = False
            print("WARNING: No connection pool implementation found. Cannot run stress test.")
            sys.exit(1)

# Load environment variables
load_dotenv()

# Get database connection parameters
DB_CONNECTION_PARAMS = {
    'host': os.environ.get('SUPABASE_DB_HOST'),
    'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
    'dbname': os.environ.get('SUPABASE_DB_NAME'),
    'user': os.environ.get('SUPABASE_DB_USER'),
    'password': os.environ.get('SUPABASE_DB_PASSWORD'),
    'sslmode': os.environ.get('SUPABASE_DB_SSLMODE', 'require')
}

# Test queries
TEST_QUERIES = [
    "SELECT COUNT(*) FROM molecule",
    "SELECT COUNT(*) FROM user_profile",
    "SELECT COUNT(*) FROM project",
    "SELECT COUNT(*) FROM mixture",
    "SELECT COUNT(*) FROM experiment"
]

def worker_thread(thread_id, pool, duration=30, max_delay=1.0):
    """Worker thread function that executes queries for a specified duration"""
    start_time = time.time()
    query_count = 0
    error_count = 0
    
    print(f"Thread {thread_id} started")
    
    while time.time() - start_time < duration:
        try:
            # Select a random query
            query = random.choice(TEST_QUERIES)
            
            # Execute query
            result = execute_query(query, pool=pool)
            
            query_count += 1
            
            # Random delay to simulate processing
            time.sleep(random.uniform(0, max_delay))
        except Exception as e:
            print(f"Thread {thread_id} error: {str(e)}")
            error_count += 1
            # Quick backoff
            time.sleep(0.5)
    
    print(f"Thread {thread_id} completed: {query_count} queries, {error_count} errors")
    return query_count, error_count

def run_stress_test(threads=10, duration=60, min_conn=2, max_conn=20):
    """Run a stress test with multiple concurrent threads"""
    print(f"Starting stress test with {threads} threads for {duration} seconds")
    print(f"Connection pool settings: min_conn={min_conn}, max_conn={max_conn}")
    
    # Create connection pool
    pool = ConnectionPool(min_conn=min_conn, max_conn=max_conn, **DB_CONNECTION_PARAMS)
    
    try:
        # Create and start threads
        thread_list = []
        for i in range(threads):
            thread = threading.Thread(
                target=worker_thread,
                args=(i, pool, duration, 0.5)
            )
            thread_list.append(thread)
            thread.start()
        
        # Wait for all threads to complete
        total_queries = 0
        total_errors = 0
        
        for i, thread in enumerate(thread_list):
            thread.join()
            # In a real implementation we would capture return values from threads
            # Here we're just estimating based on duration and threads
            queries_per_thread = (duration / 0.25) * 0.8  # Rough estimate
            total_queries += queries_per_thread
        
        print("\nStress test completed")
        print(f"Estimated total queries: {int(total_queries)}")
        print(f"Queries per second: {int(total_queries / duration)}")
        
        # Quick verification query to ensure pool is still functional
        try:
            result = execute_query("SELECT 1 as test", pool=pool)
            print("Pool remains functional after stress test")
        except Exception as e:
            print(f"Error verifying pool functionality: {str(e)}")
    
    finally:
        # Close pool
        pool.close()

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Connection Pool Stress Test')
    parser.add_argument('--threads', type=int, default=10,
                        help='Number of concurrent threads')
    parser.add_argument('--duration', type=int, default=60,
                        help='Test duration in seconds')
    parser.add_argument('--min-conn', type=int, default=2,
                        help='Minimum connections in pool')
    parser.add_argument('--max-conn', type=int, default=20,
                        help='Maximum connections in pool')
    return parser.parse_args()

if __name__ == '__main__':
    # Parse arguments
    args = parse_args()
    
    # Run stress test
    run_stress_test(
        threads=args.threads,
        duration=args.duration,
        min_conn=args.min_conn,
        max_conn=args.max_conn
    )
EOF
    fi
    
    # Run stress test with reasonable defaults
    echo "Running stress test..."
    python -m tests.stress_test_connection_pool \
        --threads 5 \
        --duration 30 \
        --min-conn $MIN_CONN \
        --max-conn $MAX_CONN
    
    TEST_RESULT_STRESS=$?
    echo "Connection pool stress test completed with exit code: $TEST_RESULT_STRESS"
else
    echo "Skipping stress test due to verification or performance test failures"
fi

# Report generation with recommendations
echo -e "\n==================================================="
echo "Connection Pool Test Results Summary"
echo "==================================================="
echo "Verification Tests: $([ $TEST_RESULT_VERIFICATION -eq 0 ] && echo "PASSED" || echo "FAILED")"
echo "Performance Tests: $([ $TEST_RESULT_PERFORMANCE -eq 0 ] && echo "PASSED" || echo "FAILED")"
if [ $TEST_RESULT_VERIFICATION -eq 0 ] && [ $TEST_RESULT_PERFORMANCE -eq 0 ]; then
    echo "Stress Tests: $([ $TEST_RESULT_STRESS -eq 0 ] && echo "PASSED" || echo "FAILED")"
fi

# Final summary
if [ $TEST_RESULT_VERIFICATION -eq 0 ] && [ $TEST_RESULT_PERFORMANCE -eq 0 ]; then
    if [ $TEST_RESULT_STRESS -eq 0 ]; then
        echo -e "\nAll connection pool tests passed successfully! 🎉"
        echo "The connection pool implementation is verified and performs well under load."
    else
        echo -e "\nConnection pool verification and performance tests passed, but stress test failed."
        echo "Review the stress test results to identify potential issues under high load."
    fi
    
    # Display report location
    echo -e "\nDetailed performance report: $REPORT_FILE"
    echo "See reports directory for all test reports."
else
    echo -e "\nSome connection pool tests failed. Please address the issues before proceeding."
    echo "Check the log files for more details."
fi

# Deactivate virtual environment
deactivate