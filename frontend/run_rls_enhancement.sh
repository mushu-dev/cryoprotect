#!/bin/bash
# Script to run the RLS policy enhancement

# Text formatting
BOLD="\033[1m"
GREEN="\033[0;32m"
RED="\033[0;31m"
YELLOW="\033[0;33m"
RESET="\033[0m"

echo -e "${BOLD}Running RLS Policy Enhancement${RESET}\n"

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

# Apply the RLS enhancements
echo -e "\n${YELLOW}Applying RLS policy enhancements...${RESET}"
python3 enhance_rls_policies.py
if [ $? -ne 0 ]; then
    echo -e "\n${RED}Failed to apply RLS policy enhancements.${RESET}"
    exit 1
else
    echo -e "${GREEN}Successfully applied RLS policy enhancements.${RESET}"
fi

# Run verification script
echo -e "\n${YELLOW}Verifying RLS policies...${RESET}"
python3 verify_rls_policies.py
if [ $? -ne 0 ]; then
    echo -e "\n${RED}RLS policy verification failed.${RESET}"
    echo -e "${YELLOW}See the generated report for details.${RESET}"
    exit 1
else
    echo -e "${GREEN}RLS policy verification successful.${RESET}"
fi

echo -e "\n${GREEN}${BOLD}RLS Policy Enhancement Complete!${RESET}"
echo -e "The RLS policies have been enhanced and verified."
echo -e "For more information, see RLS_POLICY_ENHANCEMENTS.md"
exit 0