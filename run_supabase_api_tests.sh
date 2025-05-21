#!/bin/bash

echo "Running Supabase API Integration Tests..."
echo

# Activate virtual environment if it exists
if [ -f .venv/bin/activate ]; then
    source .venv/bin/activate
fi

# Run the tests
python3 tests/run_supabase_api_tests.py

# Store the exit code
TEST_RESULT=$?

# Deactivate virtual environment if it was activated
if [ -f .venv/bin/deactivate ]; then
    deactivate
fi

# Exit with the test result
exit $TEST_RESULT