#!/bin/bash
# Script to run the service role policy consolidation

# Text formatting
BOLD="\033[1m"
GREEN="\033[0;32m"
RED="\033[0;31m"
YELLOW="\033[0;33m"
RESET="\033[0m"

echo -e "${BOLD}Running Service Role Policy Consolidation${RESET}\n"

# Check if Python is available
if ! command -v python3 &> /dev/null
then
    echo -e "${RED}Error: Python 3 is not installed or not in PATH${RESET}"
    exit 1
fi

# Check if required Python packages are installed
echo -e "${YELLOW}Checking dependencies...${RESET}"
python3 -c "import psycopg2" 2>/dev/null || { 
    echo -e "${RED}Error: psycopg2 is not installed. Run 'pip install psycopg2-binary'${RESET}"
    exit 1
}

# Check database connection
echo -e "${YELLOW}Testing database connection...${RESET}"
python3 -c "import db_utils; print('Connection successful' if db_utils.test_connection() else 'Connection failed')" || {
    echo -e "${RED}Error: Failed to connect to the database. Check your database configuration.${RESET}"
    exit 1
}

# Make script executable
chmod +x consolidate_service_role_policies.py

# Apply the service role policy consolidation
echo -e "\n${YELLOW}Applying service role policy consolidation...${RESET}"
python3 consolidate_service_role_policies.py
RESULT=$?

if [ $RESULT -ne 0 ]; then
    echo -e "\n${RED}Failed to consolidate service role policies.${RESET}"
    echo -e "${YELLOW}See the generated report for details.${RESET}"
    exit 1
else
    echo -e "\n${GREEN}${BOLD}Service Role Policy Consolidation Complete!${RESET}"
    echo -e "The service role policies have been consolidated and performance optimized."
    echo -e "Run 'python3 verify_rls_policies.py' to verify all RLS policies."
fi

exit 0