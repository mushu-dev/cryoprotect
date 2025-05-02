#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 Database Performance Testing"
echo "==============================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in the PATH."
    echo "Please install Python 3.7+ and try again."
    exit 1
fi

# Check if required packages are installed
echo "Checking required packages..."
if ! python3 -c "import supabase" &> /dev/null; then
    echo "Installing required packages..."
    pip3 install supabase psutil python-dotenv
fi

# Check if .env file exists
if [ ! -f .env ]; then
    echo "Warning: .env file not found."
    echo "Creating a template .env file. Please edit it with your Supabase credentials."
    echo "SUPABASE_URL=your_supabase_url" > .env
    echo "SUPABASE_KEY=your_supabase_service_role_key" >> .env
    echo
    echo "Please edit the .env file with your Supabase credentials and run this script again."
    exit 1
fi

echo "Starting performance tests..."
echo "This may take several minutes to complete."
echo

# Run the main test script
python3 run_all_performance_tests.py

if [ $? -ne 0 ]; then
    echo
    echo "Error: Performance tests failed."
    echo "Please check the logs for more information."
    exit 1
fi

echo
echo "==============================================================="
echo "Performance tests completed successfully!"
echo "==============================================================="
echo
echo "Reports are available in:"
echo "- comprehensive_performance_report.txt (Summary report)"
echo "- comprehensive_performance_report.json (Detailed JSON report)"
echo "- database_performance_report.txt (Database performance)"
echo "- query_plan_analysis_report.txt (Query plan analysis)"
echo "- concurrent_load_test_report.txt (Concurrent load test)"
echo
echo "See README_Performance_Testing.md for information on interpreting the results."
echo

# Ask if user wants to open the report
read -p "Would you like to open the summary report now? (y/n): " open_report

if [[ $open_report == "y" || $open_report == "Y" ]]; then
    if command -v xdg-open &> /dev/null; then
        xdg-open comprehensive_performance_report.txt  # Linux
    elif command -v open &> /dev/null; then
        open comprehensive_performance_report.txt  # macOS
    else
        echo "Could not find a program to open the report."
        echo "Please open comprehensive_performance_report.txt manually."
    fi
fi

echo
echo "Thank you for using CryoProtect v2 Database Performance Testing."
echo