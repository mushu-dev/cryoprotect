#!/bin/bash
#
# CryoProtect Production Test Suite
# This is the main driver script for production testing

set -e  # Exit on any error

# Command-line arguments
ACTION="all"
if [ $# -gt 0 ]; then
    ACTION="$1"
fi

# Script variables
TEST_RESULTS_DIR="production_test_results"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
TEST_RUN_ID="test_run_$TIMESTAMP"
TEST_OUTPUT_DIR="$TEST_RESULTS_DIR/$TEST_RUN_ID"

# Create output directory
mkdir -p "$TEST_OUTPUT_DIR"

# Log file for this test run
LOG_FILE="$TEST_OUTPUT_DIR/test_suite.log"

# Function to log messages
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# Function to display usage information
usage() {
    echo "Usage: $0 [ACTION]"
    echo "Actions:"
    echo "  setup     - Set up the container environment"
    echo "  run       - Run the tests"
    echo "  monitor   - Start the monitoring dashboard"
    echo "  cleanup   - Clean up containers and temporary files"
    echo "  all       - Perform all actions (default)"
    echo "  help      - Display this help message"
}

# Function to set up the container environment
setup_environment() {
    log "Setting up container environment..."
    
    # Make sure scripts are executable
    chmod +x production_test_container.sh run_production_tests.sh
    
    # Run setup script
    ./production_test_container.sh > "$TEST_OUTPUT_DIR/setup.log" 2>&1
    
    log "Container environment set up successfully."
}

# Function to run the tests
run_tests() {
    log "Running production tests..."
    
    # Run tests
    ./run_production_tests.sh > "$TEST_OUTPUT_DIR/tests.log" 2>&1
    
    # Copy results to output directory
    cp real_data_test_results.json "$TEST_OUTPUT_DIR/"
    cp real_data_test.log "$TEST_OUTPUT_DIR/"
    cp production_test_report.html "$TEST_OUTPUT_DIR/"
    
    log "Production tests completed. Results saved to $TEST_OUTPUT_DIR/"
}

# Function to start the monitoring dashboard
start_monitor() {
    log "Starting monitoring dashboard..."
    
    # Make monitor script executable
    chmod +x monitor_production_tests.py
    
    # Start monitor in background or foreground
    if [ "$ACTION" == "all" ]; then
        # For all action, run in background
        log "Running monitor in background mode. See $TEST_OUTPUT_DIR/monitor.log for output."
        nohup python monitor_production_tests.py > "$TEST_OUTPUT_DIR/monitor.log" 2>&1 &
        MONITOR_PID=$!
        log "Monitor started with PID $MONITOR_PID"
        # Save PID for later cleanup
        echo $MONITOR_PID > "$TEST_OUTPUT_DIR/monitor.pid"
    else
        # For monitor-only action, run in foreground
        log "Running monitor in foreground mode. Press Ctrl+C to exit."
        python monitor_production_tests.py
    fi
}

# Function to clean up
cleanup() {
    log "Cleaning up..."
    
    # Stop monitor if running
    if [ -f "$TEST_OUTPUT_DIR/monitor.pid" ]; then
        MONITOR_PID=$(cat "$TEST_OUTPUT_DIR/monitor.pid")
        if kill -0 $MONITOR_PID 2>/dev/null; then
            log "Stopping monitor process (PID $MONITOR_PID)..."
            kill $MONITOR_PID
        fi
    fi
    
    # Stop and remove containers
    log "Stopping containers..."
    podman stop cryoprotect-app cryoprotect-rdkit cryoprotect-test || true
    podman rm cryoprotect-app cryoprotect-rdkit cryoprotect-test || true
    
    log "Cleanup completed."
}

# Function to archive results
archive_results() {
    log "Archiving test results..."
    
    # Create a timestamped zip file
    ARCHIVE_FILE="production_test_results_$TIMESTAMP.zip"
    zip -r "$ARCHIVE_FILE" "$TEST_OUTPUT_DIR" > /dev/null
    
    log "Test results archived to $ARCHIVE_FILE"
}

# Main execution
case "$ACTION" in
    setup)
        log "Starting setup action..."
        setup_environment
        log "Setup action completed."
        ;;
    run)
        log "Starting run action..."
        run_tests
        log "Run action completed."
        ;;
    monitor)
        log "Starting monitor action..."
        start_monitor
        log "Monitor action completed."
        ;;
    cleanup)
        log "Starting cleanup action..."
        cleanup
        log "Cleanup action completed."
        ;;
    all)
        log "Starting full test suite..."
        setup_environment
        start_monitor
        run_tests
        cleanup
        archive_results
        log "Full test suite completed."
        ;;
    help)
        usage
        ;;
    *)
        echo "Unknown action: $ACTION"
        usage
        exit 1
        ;;
esac

log "Production test suite execution completed."