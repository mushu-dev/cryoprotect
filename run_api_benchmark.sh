#!/bin/bash

echo "Running API Performance Benchmark..."

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is not installed. Please install Python 3 and try again."
    exit 1
fi

# Check if psutil is installed
if ! python3 -c "import psutil" &> /dev/null; then
    echo "Installing psutil..."
    pip3 install psutil
fi

# Make the benchmark script executable
chmod +x benchmark_api_endpoints.py

# Run the benchmark script
python3 benchmark_api_endpoints.py "$@"

if [ $? -ne 0 ]; then
    echo "Benchmark failed with error code $?"
    exit 1
fi

echo "Benchmark completed successfully."
echo "Report saved to API_Performance_Benchmark_Report.md"