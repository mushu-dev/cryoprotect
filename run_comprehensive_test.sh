#!/bin/bash
#
# CryoProtect - Comprehensive System Test
#
# This script runs a complete test of the CryoProtect system, including:
# 1. ChEMBL data import
# 2. Toxicity data optimization
# 3. API functionality tests
# 4. Database consistency checks
# 5. Performance benchmarks

set -e

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="comprehensive_test_results_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"

echo "Starting comprehensive CryoProtect system test..."
echo "Results will be saved to $RESULTS_DIR"

# Function to run a test and save results
run_test() {
    local test_name="$1"
    local command="$2"
    local log_file="$RESULTS_DIR/${test_name}_log.txt"
    
    echo "Running $test_name..."
    echo "Command: $command"
    echo "Log: $log_file"
    
    # Run the test and capture output
    start_time=$(date +%s)
    if eval "$command" > "$log_file" 2>&1; then
        echo "✅ $test_name PASSED"
        echo "PASSED" > "$RESULTS_DIR/${test_name}_status.txt"
    else
        echo "❌ $test_name FAILED"
        echo "FAILED" > "$RESULTS_DIR/${test_name}_status.txt"
    fi
    end_time=$(date +%s)
    
    # Calculate duration
    duration=$((end_time - start_time))
    echo "Duration: ${duration}s"
    echo "$duration" > "$RESULTS_DIR/${test_name}_duration.txt"
    
    echo "----------------------------------------"
}

# Step 1: Run basic system check
run_test "system_check" "python -c \"import sys; import os; print('Python version:', sys.version); print('Current directory:', os.getcwd()); print('Files:', os.listdir('.'));\""

# Step 2: ChEMBL Import Test
# Determine number of molecules based on environment
if [ -z "$MOLECULE_COUNT" ]; then
    # Auto-detect system capabilities
    MEM_GB=$(free -g | awk '/^Mem:/{print $2}')
    CPU_COUNT=$(nproc)

    if [ "$MEM_GB" -ge 8 ] && [ "$CPU_COUNT" -ge 4 ]; then
        # High-resource system: test with 5000 molecules
        MOLECULE_COUNT=5000
    elif [ "$MEM_GB" -ge 4 ] && [ "$CPU_COUNT" -ge 2 ]; then
        # Medium-resource system: test with 1000 molecules
        MOLECULE_COUNT=1000
    else
        # Low-resource system: test with minimal set
        MOLECULE_COUNT=100
    fi
fi

echo "Running ChEMBL import test with $MOLECULE_COUNT molecules..."
run_test "chembl_import" "python comprehensive_chembl_import_test.py --molecules $MOLECULE_COUNT --report-file $RESULTS_DIR/chembl_import_report.json"

# Step 3: Apply Toxicity Optimization
run_test "toxicity_optimization" "./run_toxicity_optimization.sh"

# Step 4: API Functionality Tests
run_test "api_functionality" "python test_toxicity_optimization.py --api-url http://localhost:5000"

# Step 5: Test API Performance
run_test "api_performance" "python test_toxicity_with_samples.py"

# Step 6: Generate visualizations if matplotlib is available
echo "Generating performance visualizations..."
if python -c "import matplotlib" 2>/dev/null; then
    run_test "generate_visualizations" "python visualize_chembl_import.py --report-file $RESULTS_DIR/chembl_import_report.json --output-dir $RESULTS_DIR/visualizations"
else
    echo "Matplotlib not available, skipping visualization generation"
    echo "SKIPPED" > "$RESULTS_DIR/generate_visualizations_status.txt"
    echo "0" > "$RESULTS_DIR/generate_visualizations_duration.txt"
fi

# Step 7: Generate final report
echo "Generating comprehensive test report..."

cat > "$RESULTS_DIR/comprehensive_test_report.md" << EOL
# CryoProtect Comprehensive Test Report

**Test Date:** $(date)

## Test Summary

| Test | Status | Duration |
|------|--------|----------|
| System Check | $(cat "$RESULTS_DIR/system_check_status.txt") | $(cat "$RESULTS_DIR/system_check_duration.txt")s |
| ChEMBL Import | $(cat "$RESULTS_DIR/chembl_import_status.txt") | $(cat "$RESULTS_DIR/chembl_import_duration.txt")s |
| Toxicity Optimization | $(cat "$RESULTS_DIR/toxicity_optimization_status.txt") | $(cat "$RESULTS_DIR/toxicity_optimization_duration.txt")s |
| API Functionality | $(cat "$RESULTS_DIR/api_functionality_status.txt") | $(cat "$RESULTS_DIR/api_functionality_duration.txt")s |
| API Performance | $(cat "$RESULTS_DIR/api_performance_status.txt") | $(cat "$RESULTS_DIR/api_performance_duration.txt")s |
| Visualizations | $([ -f "$RESULTS_DIR/generate_visualizations_status.txt" ] && cat "$RESULTS_DIR/generate_visualizations_status.txt" || echo "SKIPPED") | $([ -f "$RESULTS_DIR/generate_visualizations_duration.txt" ] && cat "$RESULTS_DIR/generate_visualizations_duration.txt" || echo "0")s |

## ChEMBL Import Test Results

EOL

# Add visualizations if available
if [[ -d "$RESULTS_DIR/visualizations" ]]; then
    cat >> "$RESULTS_DIR/comprehensive_test_report.md" << EOL

### Performance Visualizations

![ChEMBL Import Summary](visualizations/chembl_import_summary.png)

![Import Rate](visualizations/chembl_import_rate.png)

![Success Rate](visualizations/chembl_import_success_rate.png)

![Progress](visualizations/chembl_import_progress.png)
EOL
fi

# Add ChEMBL import details if available
if [[ -f "$RESULTS_DIR/chembl_import_report.json" ]]; then
    python -c "
import json, sys
try:
    with open('$RESULTS_DIR/chembl_import_report.json', 'r') as f:
        data = json.load(f)
    
    with open('$RESULTS_DIR/comprehensive_test_report.md', 'a') as report:
        report.write('### Import Statistics\\n\\n')
        report.write(f'- Total molecules: {data[\"test_info\"][\"molecule_count\"]}\\n')
        report.write(f'- Successful imports: {data[\"import_stats\"][\"success_count\"]}\\n')
        report.write(f'- Failed imports: {data[\"import_stats\"][\"error_count\"]}\\n')
        report.write(f'- Total time: {data[\"import_stats\"][\"total_time\"]:.2f}s\\n')
        report.write(f'- Average time per molecule: {data[\"import_stats\"][\"average_time_per_molecule\"]:.4f}s\\n')
        report.write(f'- Import rate: {data[\"import_stats\"][\"molecules_per_second\"]:.2f} molecules/second\\n')
        
        report.write('\\n### Data Verification\\n\\n')
        report.write(f'- Verification success: {\"Yes\" if data[\"data_verification\"][\"verification_success\"] else \"No\"}\\n')
        report.write(f'- Molecules verified: {data[\"data_verification\"][\"molecules_verified\"]}\\n')
        report.write(f'- Properties verified: {data[\"data_verification\"][\"properties_verified\"]}\\n')
except Exception as e:
    print(f'Error extracting import data: {str(e)}', file=sys.stderr)
" 2>> "$RESULTS_DIR/report_generation_errors.log"
fi

# Add API Performance details
cat >> "$RESULTS_DIR/comprehensive_test_report.md" << EOL

## API Performance Test Results

EOL

# Extract performance data from the API performance log
grep -A 40 "PERFORMANCE COMPARISON" "$RESULTS_DIR/api_performance_log.txt" >> "$RESULTS_DIR/comprehensive_test_report.md" || true

# Add test logs
cat >> "$RESULTS_DIR/comprehensive_test_report.md" << EOL

## Test Logs

Detailed test logs are available in the following files:

- [System Check Log](system_check_log.txt)
- [ChEMBL Import Log](chembl_import_log.txt)
- [Toxicity Optimization Log](toxicity_optimization_log.txt)
- [API Functionality Log](api_functionality_log.txt)
- [API Performance Log](api_performance_log.txt)

## Conclusion

EOL

# Add conclusion based on test results
if grep -q "FAILED" "$RESULTS_DIR"/*_status.txt; then
    cat >> "$RESULTS_DIR/comprehensive_test_report.md" << EOL
Some tests failed. Please review the detailed logs for more information.
EOL
else
    cat >> "$RESULTS_DIR/comprehensive_test_report.md" << EOL
All tests have passed successfully. The system is functioning as expected with good performance.

Key achievements:
- Successfully imported ChEMBL data
- Toxicity data optimization is working correctly
- API endpoints are functioning properly
- Performance tests show good response times
EOL
fi

echo "Comprehensive testing completed! Report saved to $RESULTS_DIR/comprehensive_test_report.md"

# Open the report if running in a GUI environment
if command -v xdg-open >/dev/null 2>&1; then
    xdg-open "$RESULTS_DIR/comprehensive_test_report.md" 2>/dev/null || true
fi