#!/bin/bash
#
# CryoProtect - Run Full Data Pipeline Test
#
# This script tests the complete data pipeline:
# 1. From ChEMBL to Supabase database
# 2. Through API endpoints for toxicity data
# 3. With performance benchmarking

set -e

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="pipeline_test_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR/visualizations"

# Target number of molecules to import
if [ -z "$MOLECULE_COUNT" ]; then
    MOLECULE_COUNT=100
fi

echo "========================================================"
echo "  STARTING FULL DATA PIPELINE TEST"
echo "  Target: $MOLECULE_COUNT molecules"
echo "  Results will be saved to: $RESULTS_DIR"
echo "========================================================"

# Step 1: Check database connection
echo "Checking database connection..."
DB_STATUS=$(python -c "
try:
    import os
    from database.db import init_connection_pool, test_connection

    # Try to initialize connection pool
    if 'DB_HOST' in os.environ and 'DB_USER' in os.environ and 'DB_PASSWORD' in os.environ:
        init_connection_pool()
        success, message = test_connection()
        if success:
            print('CONNECTED')
        else:
            print(f'ERROR: {message}')
    else:
        print('ERROR: Missing required environment variables (DB_HOST, DB_USER, DB_PASSWORD)')
except Exception as e:
    print(f'ERROR: {str(e)}')
")

if [[ "$DB_STATUS" == CONNECTED* ]]; then
    echo "✅ Database connection available"
    HAS_DB_CONNECTION=true
else
    echo "⚠️ Database connection not available: $DB_STATUS"
    echo "Will run in test mode without database writes"
    HAS_DB_CONNECTION=false
fi

# Step 2: Check system resources
echo "Checking system resources..."
MEM_GB=$(free -g | awk '/^Mem:/{print $2}')
CPU_COUNT=$(nproc)
DISK_FREE=$(df -h . | awk 'NR==2 {print $4}')

echo "System resources:"
echo "- Memory: ${MEM_GB}GB"
echo "- CPU cores: ${CPU_COUNT}"
echo "- Free disk space: ${DISK_FREE}"
echo

# Step 3: Make script executable
chmod +x test_chembl_to_supabase_pipeline.py

# Step 4: Run the ChEMBL to Supabase pipeline test
echo "Running ChEMBL to Supabase pipeline test..."

if $HAS_DB_CONNECTION; then
    # Run with actual database connection
    ./test_chembl_to_supabase_pipeline.py \
        --molecules $MOLECULE_COUNT \
        --report-file "$RESULTS_DIR/pipeline_report.json" \
        > "$RESULTS_DIR/pipeline_log.txt" 2>&1
else
    # Run in test mode
    ./test_chembl_to_supabase_pipeline.py \
        --molecules $MOLECULE_COUNT \
        --report-file "$RESULTS_DIR/pipeline_report.json" \
        --test-mode \
        > "$RESULTS_DIR/pipeline_log.txt" 2>&1
fi

PIPELINE_EXIT_CODE=$?
echo "Pipeline test exit code: $PIPELINE_EXIT_CODE"

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo "✅ ChEMBL to Supabase pipeline test completed successfully"
else
    echo "❌ ChEMBL to Supabase pipeline test failed"
fi

# Step 5: Test API access to the inserted data
echo "Testing API access to the inserted data..."

# Check if API is already running
if pgrep -f "python app.py" > /dev/null; then
    echo "API is already running, will use existing server"
    EXISTING_API=true
else
    echo "Starting API temporarily for testing..."
    python app.py &
    API_PID=$!
    sleep 3  # Wait for the API to start
    EXISTING_API=false
fi

# Get a sample molecule ID
if $HAS_DB_CONNECTION; then
    SAMPLE_ID=$(python -c "
try:
    from database.db import execute_query
    result = execute_query('SELECT id FROM molecules ORDER BY created_at DESC LIMIT 1')
    if result and result[0]:
        print(result[0]['id'])
    else:
        import uuid
        print(uuid.uuid4())
except:
    import uuid
    print(uuid.uuid4())
")
else
    SAMPLE_ID=$(python -c "import uuid; print(uuid.uuid4())")
fi

echo "Testing API endpoints with sample ID: $SAMPLE_ID"

# Test molecule endpoint
if command -v curl &> /dev/null; then
    # Test molecule endpoint
    echo "Testing molecule API endpoint..."
    curl -s "http://localhost:5000/api/molecules/$SAMPLE_ID" -o "$RESULTS_DIR/molecule_response.json"
    
    # Test properties endpoint
    echo "Testing properties API endpoint..."
    curl -s "http://localhost:5000/api/properties/molecule/$SAMPLE_ID" -o "$RESULTS_DIR/properties_response.json"
    
    # Test toxicity endpoint
    echo "Testing toxicity API endpoint..."
    curl -s "http://localhost:5000/api/toxicity/molecule/$SAMPLE_ID" -o "$RESULTS_DIR/toxicity_response.json"
else
    # Use Python if curl is not available
    echo "Using Python to test API endpoints (curl not available)..."
    python -c "
import requests, json, sys

try:
    # Test molecule endpoint
    response = requests.get(f'http://localhost:5000/api/molecules/$SAMPLE_ID')
    with open('$RESULTS_DIR/molecule_response.json', 'w') as f:
        json.dump(response.json() if response.status_code == 200 else {'error': response.text}, f, indent=2)
    
    # Test properties endpoint
    response = requests.get(f'http://localhost:5000/api/properties/molecule/$SAMPLE_ID')
    with open('$RESULTS_DIR/properties_response.json', 'w') as f:
        json.dump(response.json() if response.status_code == 200 else {'error': response.text}, f, indent=2)
    
    # Test toxicity endpoint
    response = requests.get(f'http://localhost:5000/api/toxicity/molecule/$SAMPLE_ID')
    with open('$RESULTS_DIR/toxicity_response.json', 'w') as f:
        json.dump(response.json() if response.status_code == 200 else {'error': response.text}, f, indent=2)
        
    print('API tests completed. Responses saved to files.')
except Exception as e:
    print(f'Error testing API: {str(e)}', file=sys.stderr)
"
fi

# Clean up if we started the API
if [[ "$EXISTING_API" == "false" ]]; then
    echo "Stopping temporary API server..."
    kill $API_PID
fi

# Step 6: Apply toxicity optimization and test performance
echo "Applying toxicity optimization and testing performance..."
bash run_toxicity_optimization.sh > "$RESULTS_DIR/toxicity_optimization_log.txt" 2>&1
TOXICITY_EXIT_CODE=$?

if [ $TOXICITY_EXIT_CODE -eq 0 ]; then
    echo "✅ Toxicity optimization applied successfully"
else
    echo "⚠️ Toxicity optimization had some issues, see log for details"
fi

# Step 7: Test toxicity API performance
echo "Testing toxicity API performance..."
python test_toxicity_with_samples.py > "$RESULTS_DIR/toxicity_performance_log.txt" 2>&1
PERFORMANCE_EXIT_CODE=$?

if [ $PERFORMANCE_EXIT_CODE -eq 0 ]; then
    echo "✅ Toxicity API performance test completed successfully"
else
    echo "⚠️ Toxicity API performance test had some issues, see log for details"
fi

# Step 8: Generate visualizations if available
if [ -f "$RESULTS_DIR/pipeline_report.json" ]; then
    echo "Generating pipeline visualizations..."
    
    # Install matplotlib if needed and available
    if ! python -c "import matplotlib" 2>/dev/null; then
        echo "Matplotlib not available, attempting to install..."
        pip install matplotlib
    fi
    
    # Generate visualizations
    if python -c "import matplotlib" 2>/dev/null; then
        python -c "
import json, os, sys
import matplotlib.pyplot as plt

try:
    # Load data
    with open('$RESULTS_DIR/pipeline_report.json', 'r') as f:
        data = json.load(f)
    
    # Create output directory
    os.makedirs('$RESULTS_DIR/visualizations', exist_ok=True)
    
    # Plot 1: Pipeline Time Breakdown
    plt.figure(figsize=(10, 6))
    labels = ['ChEMBL Fetch', 'Supabase Insert', 'Other']
    times = [
        data['pipeline_stats']['chembl_fetch_time'] * data['pipeline_stats']['success_count'],
        data['pipeline_stats']['supabase_insert_time'] * data['pipeline_stats']['success_count'],
        data['pipeline_stats']['total_time'] - 
        (data['pipeline_stats']['chembl_fetch_time'] + data['pipeline_stats']['supabase_insert_time']) * data['pipeline_stats']['success_count']
    ]
    plt.pie(times, labels=labels, autopct='%1.1f%%')
    plt.title('Pipeline Time Breakdown')
    plt.savefig('$RESULTS_DIR/visualizations/pipeline_time_breakdown.png')
    
    # Plot 2: Success Rate
    plt.figure(figsize=(10, 6))
    labels = ['Success', 'Error']
    counts = [data['pipeline_stats']['success_count'], data['pipeline_stats']['error_count']]
    plt.pie(counts, labels=labels, colors=['#4CAF50', '#F44336'], autopct='%1.1f%%')
    plt.title('Pipeline Success Rate')
    plt.savefig('$RESULTS_DIR/visualizations/pipeline_success_rate.png')
    
    # Plot 3: Database Changes
    if 'database_stats' in data and 'difference' in data['database_stats'] and 'table_counts' in data['database_stats']['difference']:
        diff = data['database_stats']['difference']['table_counts']
        plt.figure(figsize=(10, 6))
        tables = list(diff.keys())
        counts = [diff[table] for table in tables]
        plt.bar(tables, counts, color='#2196F3')
        plt.title('Database Changes')
        plt.ylabel('Rows Added')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig('$RESULTS_DIR/visualizations/database_changes.png')
    
    print('Visualizations generated successfully')
except Exception as e:
    print(f'Error generating visualizations: {str(e)}', file=sys.stderr)
"
    else
        echo "Skipping visualizations due to missing matplotlib"
    fi
fi

# Step 9: Generate summary report
echo "Generating summary report..."
cat > "$RESULTS_DIR/full_pipeline_test_report.md" << EOL
# CryoProtect Full Data Pipeline Test Report

**Test Date:** $(date)
**Target Molecules:** ${MOLECULE_COUNT}

## System Information

- Memory: ${MEM_GB}GB
- CPU Cores: ${CPU_COUNT}
- Free Disk Space: ${DISK_FREE}
- Database Connection: $([ "$HAS_DB_CONNECTION" == "true" ] && echo "Available" || echo "Not available")

## Test Results

| Component | Status | Notes |
|-----------|--------|-------|
| ChEMBL to Supabase Pipeline | $([ $PIPELINE_EXIT_CODE -eq 0 ] && echo "✅ PASSED" || echo "❌ FAILED") | See [pipeline_report.md](pipeline_report.md) for details |
| Toxicity Optimization | $([ $TOXICITY_EXIT_CODE -eq 0 ] && echo "✅ PASSED" || echo "⚠️ PARTIAL") | See [toxicity_optimization_log.txt](toxicity_optimization_log.txt) for details |
| Toxicity API Performance | $([ $PERFORMANCE_EXIT_CODE -eq 0 ] && echo "✅ PASSED" || echo "⚠️ PARTIAL") | See [toxicity_performance_log.txt](toxicity_performance_log.txt) for details |

## Test Details

EOL

# Add pipeline details if available
if [ -f "$RESULTS_DIR/pipeline_report.json" ]; then
python -c "
import json, sys
try:
    with open('$RESULTS_DIR/pipeline_report.json', 'r') as f:
        data = json.load(f)
    
    with open('$RESULTS_DIR/full_pipeline_test_report.md', 'a') as report:
        # Pipeline stats
        stats = data['pipeline_stats']
        report.write('### ChEMBL to Supabase Pipeline\\n\\n')
        report.write(f'- Total molecules: {data[\"test_info\"][\"molecule_count\"]}\\n')
        report.write(f'- Success rate: {stats[\"success_count\"]}/{data[\"test_info\"][\"molecule_count\"]} ({stats[\"success_count\"]/data[\"test_info\"][\"molecule_count\"]*100:.1f}%)\\n')
        report.write(f'- Total time: {stats[\"total_time\"]:.2f} seconds\\n')
        report.write(f'- Overall rate: {stats[\"molecules_per_second\"]:.2f} molecules/second\\n')
        report.write(f'- Average ChEMBL fetch time: {stats[\"chembl_fetch_time\"]*1000:.2f} ms per molecule\\n')
        report.write(f'- Average Supabase insert time: {stats[\"supabase_insert_time\"]*1000:.2f} ms per molecule\\n\\n')
        
        # Database changes
        if 'database_stats' in data and 'difference' in data['database_stats'] and 'table_counts' in data['database_stats']['difference']:
            diff = data['database_stats']['difference']['table_counts']
            report.write('#### Database Changes\\n\\n')
            for table, count in diff.items():
                report.write(f'- {table}: +{count} rows\\n')
            report.write('\\n')
        
        # Add visualizations if available
        if '$RESULTS_DIR/visualizations/pipeline_time_breakdown.png' and '$RESULTS_DIR/visualizations/pipeline_success_rate.png':
            report.write('#### Pipeline Visualizations\\n\\n')
            report.write('![Pipeline Time Breakdown](visualizations/pipeline_time_breakdown.png)\\n\\n')
            report.write('![Pipeline Success Rate](visualizations/pipeline_success_rate.png)\\n\\n')
            if '$RESULTS_DIR/visualizations/database_changes.png':
                report.write('![Database Changes](visualizations/database_changes.png)\\n\\n')
except Exception as e:
    print(f'Error extracting pipeline data: {str(e)}', file=sys.stderr)
" 2>> "$RESULTS_DIR/report_generation_errors.log"
fi

# Add toxicity performance details if available
if grep -q "PERFORMANCE COMPARISON" "$RESULTS_DIR/toxicity_performance_log.txt"; then
    echo "Adding toxicity performance metrics to report..."
    cat >> "$RESULTS_DIR/full_pipeline_test_report.md" << EOL

### Toxicity API Performance

EOL
    grep -A 15 "PERFORMANCE COMPARISON" "$RESULTS_DIR/toxicity_performance_log.txt" | sed 's/=/\-/g' >> "$RESULTS_DIR/full_pipeline_test_report.md"
fi

# Add conclusion
cat >> "$RESULTS_DIR/full_pipeline_test_report.md" << EOL

## Conclusion

$([ $PIPELINE_EXIT_CODE -eq 0 ] && [ $TOXICITY_EXIT_CODE -eq 0 ] && [ $PERFORMANCE_EXIT_CODE -eq 0 ] && echo "✅ All pipeline tests passed successfully. The system is functioning properly." || echo "⚠️ Some tests had issues. Please review the detailed logs for more information.")

## Next Steps

1. Review the detailed reports for any issues
2. Deploy the optimization to production if all tests passed
3. Set up monitoring for the database and API endpoints
4. Update documentation with the optimization details
EOL

# Step 10: Open the report if running in a GUI environment
if command -v xdg-open >/dev/null 2>&1; then
    echo "Opening summary report..."
    xdg-open "$RESULTS_DIR/full_pipeline_test_report.md" 2>/dev/null || true
fi

echo "========================================================"
echo "  FULL DATA PIPELINE TEST COMPLETED"
echo "  Results saved to: $RESULTS_DIR"
echo "  Summary report: $RESULTS_DIR/full_pipeline_test_report.md"
echo "========================================================"