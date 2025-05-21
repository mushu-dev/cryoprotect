#!/bin/bash
# Environment Variable Validation Script for CryoProtect v2
#
# This script validates that all required environment variables are set by:
# 1. Parsing the .env.template file to extract variable names
# 2. Checking if each variable is set in the current environment
# 3. Reporting any missing variables
# 4. Exiting with appropriate status code

# Default values
TEMPLATE_FILE=".env.template"
CHECK_ALL=false
VERBOSE=false
SECTION_FILTER=""

# Function to display usage information
usage() {
  echo "Usage: $0 [options]"
  echo ""
  echo "Options:"
  echo "  -t, --template PATH     Path to the .env.template file (default: .env.template)"
  echo "  -s, --section SECTION   Only check variables in the specified section"
  echo "  -a, --all               Check all variables, not just required ones"
  echo "  -v, --verbose           Show more detailed output"
  echo "  -h, --help              Show this help message and exit"
  exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -t|--template)
      TEMPLATE_FILE="$2"
      shift 2
      ;;
    -s|--section)
      SECTION_FILTER="$2"
      shift 2
      ;;
    -a|--all)
      CHECK_ALL=true
      shift
      ;;
    -v|--verbose)
      VERBOSE=true
      shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

# Check if template file exists
if [ ! -f "$TEMPLATE_FILE" ]; then
  echo "Error: Template file '$TEMPLATE_FILE' not found."
  exit 1
fi

echo "Validating environment variables for CryoProtect v2 using template: $TEMPLATE_FILE"
if [ -n "$SECTION_FILTER" ]; then
  echo "Checking only variables in section: $SECTION_FILTER"
fi
if [ "$CHECK_ALL" = true ]; then
  echo "Checking all variables (required and optional)"
fi

# Initialize counters and arrays
TOTAL_VARS=0
REQUIRED_VARS=0
OPTIONAL_VARS=0
MISSING_REQUIRED=()
MISSING_OPTIONAL=()
CURRENT_SECTION="GENERAL"

# Process the template file
while IFS= read -r line; do
  # Trim whitespace
  line=$(echo "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  
  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi
  
  # Check if this is a section header
  if [[ "$line" =~ ^#[[:space:]]*=====[[:space:]]*(.*)[[:space:]]*=====[[:space:]]*$ ]]; then
    CURRENT_SECTION=$(echo "${BASH_REMATCH[1]}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    continue
  fi
  
  # Skip comment-only lines
  if [[ "$line" =~ ^#.*$ ]]; then
    continue
  fi
  
  # Extract variable name from lines with equals sign
  if [[ "$line" =~ ^([A-Za-z0-9_]+)= ]]; then
    VAR_NAME="${BASH_REMATCH[1]}"
    TOTAL_VARS=$((TOTAL_VARS + 1))
    
    # Skip if we're filtering by section and this isn't in the target section
    if [ -n "$SECTION_FILTER" ] && [ "$SECTION_FILTER" != "$CURRENT_SECTION" ]; then
      continue
    fi
    
    # Check if this variable is required by looking at nearby comments
    IS_REQUIRED=false
    
    # Look at the previous 5 lines for comments about this variable
    prev_lines=$(grep -B 5 "^$VAR_NAME=" "$TEMPLATE_FILE" | grep "^#")
    if echo "$prev_lines" | grep -i "required" > /dev/null && ! echo "$prev_lines" | grep -i "optional" > /dev/null; then
      IS_REQUIRED=true
      REQUIRED_VARS=$((REQUIRED_VARS + 1))
    else
      # Apply heuristics for certain variables
      if [[ "$CURRENT_SECTION" == "DATABASE CONFIGURATION" && ("$VAR_NAME" == "SUPABASE_URL" || "$VAR_NAME" == "SUPABASE_KEY") ]]; then
        IS_REQUIRED=true
        REQUIRED_VARS=$((REQUIRED_VARS + 1))
      elif [[ "$CURRENT_SECTION" == "APPLICATION CONFIGURATION" && "$VAR_NAME" == "SECRET_KEY" ]]; then
        IS_REQUIRED=true
        REQUIRED_VARS=$((REQUIRED_VARS + 1))
      else
        OPTIONAL_VARS=$((OPTIONAL_VARS + 1))
      fi
    fi
    
    # Check if the variable is set in the environment
    if [ -z "${!VAR_NAME}" ]; then
      if [ "$IS_REQUIRED" = true ]; then
        MISSING_REQUIRED+=("$VAR_NAME:$CURRENT_SECTION")
      elif [ "$CHECK_ALL" = true ]; then
        MISSING_OPTIONAL+=("$VAR_NAME:$CURRENT_SECTION")
      fi
    fi
  fi
done < "$TEMPLATE_FILE"

if [ "$VERBOSE" = true ]; then
  echo ""
  echo "Found $TOTAL_VARS variables in template file"
  echo "  - $REQUIRED_VARS required variables"
  echo "  - $OPTIONAL_VARS optional variables"
fi

# Report results
EXIT_CODE=0

if [ ${#MISSING_REQUIRED[@]} -gt 0 ]; then
  echo ""
  echo "ERROR: The following required environment variables are missing:"
  
  # Group by section
  declare -A SECTIONS
  for var in "${MISSING_REQUIRED[@]}"; do
    IFS=':' read -r name section <<< "$var"
    SECTIONS["$section"]+="$name "
  done
  
  # Print missing variables by section
  for section in "${!SECTIONS[@]}"; do
    echo ""
    echo "$section:"
    for var in ${SECTIONS["$section"]}; do
      echo "  - $var"
    done
  done
  
  EXIT_CODE=1
fi

if [ "$CHECK_ALL" = true ] && [ ${#MISSING_OPTIONAL[@]} -gt 0 ]; then
  echo ""
  echo "WARNING: The following optional environment variables are not set:"
  
  # Group by section
  declare -A SECTIONS
  for var in "${MISSING_OPTIONAL[@]}"; do
    IFS=':' read -r name section <<< "$var"
    SECTIONS["$section"]+="$name "
  done
  
  # Print missing variables by section
  for section in "${!SECTIONS[@]}"; do
    echo ""
    echo "$section:"
    for var in ${SECTIONS["$section"]}; do
      echo "  - $var"
    done
  done
fi

if [ $EXIT_CODE -eq 1 ]; then
  echo ""
  echo "Please set these required variables in your environment or .env file before running the application."
else
  echo "All required environment variables are set."
fi

exit $EXIT_CODE