#!/bin/bash
#
# CryoProtect - Run Large Scale ChEMBL Import
#
# This script performs a full large-scale import of ChEMBL molecules
# and generates comprehensive performance reports.

set -e

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="chembl_large_import_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR/visualizations"
mkdir -p "$RESULTS_DIR/checkpoints"

# Target number of molecules to import
if [ -z "$MOLECULE_COUNT" ]; then
    MOLECULE_COUNT=10000
fi

echo "========================================================"
echo "  STARTING LARGE-SCALE CHEMBL IMPORT TEST"
echo "  Target: $MOLECULE_COUNT molecules"
echo "  Results will be saved to: $RESULTS_DIR"
echo "========================================================"

# Step 1: Check system resources
echo "Checking system resources..."
MEM_GB=$(free -g | awk '/^Mem:/{print $2}')
CPU_COUNT=$(nproc)
DISK_FREE=$(df -h . | awk 'NR==2 {print $4}')

echo "System resources:"
echo "- Memory: ${MEM_GB}GB"
echo "- CPU cores: ${CPU_COUNT}"
echo "- Free disk space: ${DISK_FREE}"
echo

# Step 2: Create cache directory if it doesn't exist
CACHE_DIR="chembl_cache_large"
mkdir -p "$CACHE_DIR"
echo "Using cache directory: $CACHE_DIR"

# Step 3: Run the import with progress monitoring
echo "Starting import of $MOLECULE_COUNT molecules..."
echo "This may take some time, please be patient."
echo

# Start time
START_TIME=$(date +%s)

# Run the import with checkpoint files going to results directory
python comprehensive_chembl_import_test.py \
    --molecules $MOLECULE_COUNT \
    --cache-dir "$CACHE_DIR" \
    --report-file "$RESULTS_DIR/chembl_import_report.json" \
    --checkpoint-dir "$RESULTS_DIR/checkpoints" \
    > "$RESULTS_DIR/import_log.txt" 2>&1 &

IMPORT_PID=$!

# Monitor progress
while kill -0 $IMPORT_PID 2>/dev/null; do
    if [ -f "$RESULTS_DIR/checkpoints/progress.json" ]; then
        # Extract progress information
        PROGRESS=$(grep -o '"progress_percentage":[0-9.]*' "$RESULTS_DIR/checkpoints/progress.json" | cut -d':' -f2)
        ELAPSED=$(grep -o '"elapsed_time":[0-9.]*' "$RESULTS_DIR/checkpoints/progress.json" | cut -d':' -f2)
        SUCCESS=$(grep -o '"success_count":[0-9]*' "$RESULTS_DIR/checkpoints/progress.json" | cut -d':' -f2)
        ERROR=$(grep -o '"error_count":[0-9]*' "$RESULTS_DIR/checkpoints/progress.json" | cut -d':' -f2)
        RATE=$(grep -o '"molecules_per_second":[0-9.]*' "$RESULTS_DIR/checkpoints/progress.json" | cut -d':' -f2)
        
        # Print progress
        echo -ne "Progress: ${PROGRESS:-0}% | Elapsed: ${ELAPSED:-0}s | Success: ${SUCCESS:-0} | Errors: ${ERROR:-0} | Rate: ${RATE:-0} mol/s \r"
    else
        echo -ne "Waiting for progress information... \r"
    fi
    sleep 2
done

# Wait for the import to finish
wait $IMPORT_PID
EXIT_CODE=$?

# End time
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo
echo "Import process finished in ${TOTAL_TIME} seconds"
echo "Exit code: $EXIT_CODE"
echo

# Step 4: Generate visualizations
echo "Generating performance visualizations..."
python visualize_chembl_import.py \
    --report-file "$RESULTS_DIR/chembl_import_report.json" \
    --output-dir "$RESULTS_DIR/visualizations" \
    > "$RESULTS_DIR/visualization_log.txt" 2>&1

# Step 5: Generate summary report
echo "Generating summary report..."
cat > "$RESULTS_DIR/summary_report.md" << EOL
# ChEMBL Large-Scale Import Test Results

**Test Date:** $(date)
**Duration:** ${TOTAL_TIME} seconds
**Target Molecules:** ${MOLECULE_COUNT}

## System Information

- Memory: ${MEM_GB}GB
- CPU Cores: ${CPU_COUNT}
- Free Disk Space: ${DISK_FREE}
- Cache Directory: ${CACHE_DIR}

## Import Performance

EOL

# Extract key metrics from the report file
python -c "
import json, sys
try:
    with open('$RESULTS_DIR/chembl_import_report.json', 'r') as f:
        data = json.load(f)
    
    with open('$RESULTS_DIR/summary_report.md', 'a') as report:
        # Import stats
        stats = data['import_stats']
        report.write(f'- Total molecules: {data[\"test_info\"][\"molecule_count\"]}\\n')
        report.write(f'- Successfully imported: {stats[\"success_count\"]}\\n')
        report.write(f'- Failed imports: {stats[\"error_count\"]}\\n')
        report.write(f'- Success rate: {stats[\"success_count\"]/data[\"test_info\"][\"molecule_count\"]*100:.2f}%\\n')
        report.write(f'- Total properties: {stats[\"total_properties_count\"]}\\n')
        report.write(f'- Import time: {stats[\"total_time\"]:.2f} seconds\\n')
        report.write(f'- Average time per molecule: {stats[\"average_time_per_molecule\"]*1000:.2f} ms\\n')
        report.write(f'- Import rate: {stats[\"molecules_per_second\"]:.2f} molecules/second\\n')
        
        # Performance metrics if available
        if 'performance_summary' in data:
            summary = data['performance_summary']
            report.write('\\n## Performance Stability\\n\\n')
            report.write(f'- Average rate: {summary[\"average_rate\"]:.2f} molecules/second\\n')
            report.write(f'- Minimum rate: {summary[\"min_rate\"]:.2f} molecules/second\\n')
            report.write(f'- Maximum rate: {summary[\"max_rate\"]:.2f} molecules/second\\n')
            report.write(f'- Rate stability: {summary[\"rate_stability\"]*100:.2f}%\\n')
        
        # Analysis if available
        if 'analysis' in data:
            analysis = data['analysis']
            report.write('\\n## Scaling Analysis\\n\\n')
            report.write(f'- Time per 1000 molecules: {analysis[\"time_per_1000_molecules\"]:.2f} seconds\\n')
            report.write(f'- Estimated time for 100K molecules: {analysis[\"time_per_1000_molecules\"]*100/60:.2f} minutes\\n')
            
        # Add visualization references
        report.write('\\n## Visualizations\\n\\n')
        report.write('![Import Summary](visualizations/chembl_import_summary.png)\\n\\n')
        report.write('![Import Rate](visualizations/chembl_import_rate.png)\\n\\n')
        report.write('![Success Rate](visualizations/chembl_import_success_rate.png)\\n\\n')
        report.write('![Progress](visualizations/chembl_import_progress.png)\\n')
        
except Exception as e:
    print(f'Error extracting data: {str(e)}', file=sys.stderr)
" 2>> "$RESULTS_DIR/report_generation_errors.log"

# Step 6: Open the report if running in a GUI environment
if command -v xdg-open >/dev/null 2>&1; then
    echo "Opening summary report..."
    xdg-open "$RESULTS_DIR/summary_report.md" 2>/dev/null || true
fi

echo "========================================================"
echo "  LARGE-SCALE CHEMBL IMPORT TEST COMPLETED"
echo "  Results saved to: $RESULTS_DIR"
echo "  Summary report: $RESULTS_DIR/summary_report.md"
echo "========================================================"