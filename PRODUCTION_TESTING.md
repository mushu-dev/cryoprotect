# CryoProtect Production Testing Guide

This guide explains how to use the production testing tools for CryoProtect. These tools allow you to verify that the system works correctly with real data in a production-like environment.

## Test Components

The production testing system consists of the following components:

1. **Container Environment**
   - `cryoprotect-app`: The main application container
   - `cryoprotect-rdkit`: The RDKit service container
   - `cryoprotect-test`: A dedicated test container

2. **Test Data**
   - Core cryoprotectants dataset
   - Edge cases dataset
   - Mixtures dataset

3. **Test Scripts**
   - `production_test_suite.sh`: Main test driver
   - `production_test_container.sh`: Container setup script
   - `run_production_tests.sh`: Test execution script
   - `monitor_production_tests.py`: Real-time monitoring dashboard
   - `production_workflow_test.py`: Core test implementation

## Running the Tests

### Option 1: Full Test Suite

To run the complete test suite (setup, monitor, run, cleanup), simply execute:

```bash
./production_test_suite.sh
```

This will:
1. Set up the container environment
2. Start the monitoring dashboard
3. Run all tests
4. Clean up containers and temporary files
5. Generate an HTML report

### Option 2: Individual Steps

You can also run individual steps of the test process:

```bash
# Set up the container environment
./production_test_suite.sh setup

# Run the tests
./production_test_suite.sh run

# Start the monitoring dashboard
./production_test_suite.sh monitor

# Clean up after testing
./production_test_suite.sh cleanup
```

## Test Results

After running the tests, you'll find the results in the `production_test_results` directory:

- `production_test_results/test_run_TIMESTAMP/`: Directory containing all test artifacts
  - `test_suite.log`: Main test log
  - `setup.log`: Container setup log
  - `tests.log`: Test execution log
  - `monitor.log`: Monitor output log
  - `production_test_report.html`: HTML report with test results
  - `production_workflow_results.json`: Raw test results in JSON format

Additionally, a ZIP archive `production_test_results_TIMESTAMP.zip` is created with all test artifacts for easy sharing.

## Monitoring Dashboard

The monitoring dashboard provides real-time visibility into:

- Container resource usage (CPU, memory)
- API endpoint performance
- Test progress and status
- Most recent log events

To start the dashboard:

```bash
./monitor_production_tests.py
```

## Customizing Tests

You can customize the production tests by:

1. Modifying `production_workflow_test.py` to add or change test cases
2. Adjusting performance thresholds in the `analyze_performance()` method
3. Adding new test data in the `tests/test_data/` directory

## Performance Benchmarks

The production tests include performance benchmarks for key operations:

- Molecule property calculation
- Molecule visualization
- Mixture creation and analysis
- Search functionality (by name and structure)
- Batch processing

These benchmarks help identify performance regressions before deploying to production.

## Troubleshooting

### Container Issues

If you encounter container-related issues:

```bash
# Check container status
podman ps -a

# View container logs
podman logs cryoprotect-app
podman logs cryoprotect-rdkit
podman logs cryoprotect-test

# Restart containers
podman restart cryoprotect-app cryoprotect-rdkit cryoprotect-test
```

### Network Issues

If containers cannot communicate:

```bash
# Check network
podman network inspect cryoprotect-net

# Verify connectivity from test container
podman exec cryoprotect-test curl -s http://cryoprotect-app:5000/health
podman exec cryoprotect-test curl -s http://cryoprotect-rdkit:5000/health
```

### Test Data Issues

If test data is missing or corrupt:

```bash
# Verify test data files
ls -la tests/test_data/

# Check file contents
cat tests/test_data/core_cryoprotectants.json | head
```

## Continuous Integration

The production tests can be integrated into CI/CD pipelines using the following command:

```bash
./production_test_suite.sh all
```

The script will exit with code 0 if all tests pass, or non-zero if any test fails.