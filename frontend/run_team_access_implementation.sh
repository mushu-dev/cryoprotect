#!/bin/bash
# Script to implement team-based access control

# Text formatting
BOLD="\033[1m"
GREEN="\033[0;32m"
RED="\033[0;31m"
YELLOW="\033[0;33m"
RESET="\033[0m"

echo -e "${BOLD}Implementing Team-Based Access Control${RESET}\n"

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
chmod +x implement_team_based_access.py

# Apply the team-based access control implementation
echo -e "\n${YELLOW}Implementing team-based access control...${RESET}"
python3 implement_team_based_access.py

RESULT=$?

if [ $RESULT -ne 0 ]; then
    echo -e "\n${RED}Failed to implement team-based access control.${RESET}"
    echo -e "${YELLOW}See the generated report for details.${RESET}"
    exit 1
fi

# Create test data
echo -e "\n${YELLOW}Creating test data to validate team-based access...${RESET}"
python3 implement_team_based_access.py --create-data
RESULT=$?

if [ $RESULT -ne 0 ]; then
    echo -e "\n${RED}Failed to create test data.${RESET}"
    exit 1
fi

echo -e "\n${GREEN}${BOLD}Team-Based Access Control Implementation Complete!${RESET}"
echo -e "The team-based access control has been implemented and verified."
echo -e "For more information, see TEAM_BASED_ACCESS_CONTROL.md"
exit 0