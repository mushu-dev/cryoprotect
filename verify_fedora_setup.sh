#!/bin/bash

# CryoProtect Fedora Setup Verification Script
# This script checks if the CryoProtect environment is properly set up on Fedora

echo "======================================================="
echo "           CryoProtect Fedora Setup Verification       "
echo "======================================================="

# Check if we're running on Fedora
echo -n "Checking if running on Fedora... "
if [ -f /etc/fedora-release ]; then
    FEDORA_VERSION=$(cat /etc/fedora-release | grep -oE '[0-9]+')
    echo "Yes (Fedora $FEDORA_VERSION)"
else
    echo "No"
    echo "WARNING: This system is not running Fedora. Some checks may not be relevant."
fi

# Check if SELinux is enabled
echo -n "Checking SELinux status... "
if command -v getenforce > /dev/null; then
    SELINUX_STATUS=$(getenforce)
    echo "$SELINUX_STATUS"
    
    if [ "$SELINUX_STATUS" == "Enforcing" ]; then
        echo "  - SELinux is enforcing - make sure contexts are properly set"
    fi
else
    echo "Not available"
fi

# Check if Python is installed and version
echo -n "Checking Python installation... "
if command -v python > /dev/null; then
    PYTHON_VERSION=$(python --version 2>&1)
    echo "$PYTHON_VERSION"
else
    echo "Not found"
    echo "ERROR: Python is not installed. Please install Python 3.x."
    exit 1
fi

# Check if the .env file exists
echo -n "Checking for .env file... "
if [ -f .env ]; then
    echo "Found"
    
    # Check if SUPABASE credentials are set
    if grep -q "SUPABASE_URL" .env && grep -q "SUPABASE_KEY" .env; then
        echo "  - Supabase credentials found in .env"
    else
        echo "  - WARNING: Supabase credentials not found in .env"
    fi
else
    echo "Not found"
    echo "WARNING: .env file not found. You will need to create one with Supabase credentials."
fi

# Check for Python virtual environment
echo -n "Checking for Python virtual environment... "
if [ -d "quick_env" ]; then
    echo "Found"
    
    # Check if the virtual environment has the basic requirements
    if [ -f "quick_env/bin/activate" ]; then
        echo "  - Virtual environment appears valid"
    else
        echo "  - WARNING: Virtual environment may be incomplete"
    fi
else
    echo "Not found"
    echo "INFO: Python virtual environment not found. It will be created when you run the minimal app."
fi

# Check if simplified_app.py exists
echo -n "Checking for simplified app... "
if [ -f "simplified_app.py" ]; then
    echo "Found"
else
    echo "Not found"
    echo "WARNING: simplified_app.py not found. The minimal app test may not work."
fi

# Check connectivity to Supabase if credentials exist
echo -n "Testing Supabase connectivity... "
if [ -f .env ] && [ -f "simplified_app.py" ]; then
    # Export environment variables from .env
    export $(grep -v '^#' .env | xargs) 2>/dev/null
    
    # Check if curl is installed
    if command -v curl > /dev/null; then
        # Start the app in the background if it's not already running
        APP_RUNNING=$(ps aux | grep "python simplified_app.py" | grep -v grep)
        
        if [ -z "$APP_RUNNING" ]; then
            echo "Starting app temporarily..."
            
            # Check if we have a virtual environment
            if [ -f "quick_env/bin/activate" ]; then
                source quick_env/bin/activate
                nohup python simplified_app.py > simplified_app.log 2>&1 &
                APP_PID=$!
                
                # Wait for the app to start
                sleep 5
                
                # Check if still running
                if ps -p $APP_PID > /dev/null; then
                    echo "  - App started with PID $APP_PID"
                else
                    echo "  - Failed to start app. Check simplified_app.log"
                    exit 1
                fi
            else
                echo "  - Cannot test connectivity - virtual environment not set up"
                echo "  - Run './run_simplified_app.sh' first to set up environment"
                exit 1
            fi
        else
            echo "App already running..."
        fi
        
        # Test the connection
        RESPONSE=$(curl -s http://localhost:5000/test-supabase)
        
        # Check if response contains "success"
        if echo "$RESPONSE" | grep -q "success"; then
            echo "  - Supabase connection successful!"
        else
            echo "  - Supabase connection failed. Response:"
            echo "$RESPONSE"
        fi
        
        # If we started the app, kill it
        if [ -n "$APP_PID" ]; then
            kill $APP_PID
            echo "  - Stopped temporary app instance"
        fi
    else
        echo "Unable to test - curl not installed"
    fi
else
    echo "Skipped (missing .env or app)"
fi

# Check RDKit installation
echo -n "Checking RDKit installation... "
if [ -f "quick_env/bin/activate" ] && [ -f "check_rdkit.py" ]; then
    source quick_env/bin/activate
    # Use our dedicated check_rdkit.py script
    RDKIT_CHECK=$(python check_rdkit.py 2>&1)
    RDKIT_STATUS=$?
    deactivate

    if [ $RDKIT_STATUS -eq 0 ]; then
        echo "$RDKIT_CHECK"
    else
        echo "Not installed"
        echo "INFO: RDKit is not installed. Install with 'pip install rdkit' or 'conda install -c conda-forge rdkit'"
    fi
else
    if [ ! -f "check_rdkit.py" ]; then
        echo "Skipped (check_rdkit.py not found)"
    else
        echo "Skipped (no environment)"
    fi
fi

echo "======================================================="
echo "               Verification Summary                    "
echo "======================================================="
echo ""
echo "To complete the setup:"
echo ""
echo "1. Run the simplified app for testing:"
echo "   ./run_simplified_app.sh"
echo ""
echo "2. Verify all endpoints work at http://localhost:5000/"
echo ""
echo "3. Check DATABASE_SCHEMA_SUMMARY.md for database structure"
echo ""
echo "4. Check FEDORA_SETUP_SUMMARY.md for next steps"
echo ""
echo "Useful verification tools:"
echo "- check_rdkit.py: Verify RDKit installation"
echo "- simplified_app.py: Test basic functionality"
echo "- run_simplified_app.sh: Run the minimal test app"
echo ""
echo "======================================================="