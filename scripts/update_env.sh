#!/bin/bash
# Updates the .env file with Supabase credentials

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

ENV_FILE="../.env"

echo -e "${BLUE}CryoProtect Environment Configuration${NC}"
echo "This script will update your .env file with Supabase credentials"
echo "==========================================================================="

# Check if .env file exists
if [ ! -f "$ENV_FILE" ]; then
  echo -e "${YELLOW}No .env file found at $ENV_FILE. Creating a new one.${NC}"
  touch "$ENV_FILE"
fi

# Function to get existing value or default
get_value() {
  local key=$1
  local default=$2
  local current_value=$(grep "^$key=" "$ENV_FILE" | cut -d= -f2-)
  
  if [ -n "$current_value" ]; then
    echo "$current_value"
  else
    echo "$default"
  fi
}

# Get existing values
CURRENT_URL=$(get_value "SUPABASE_URL" "https://your-project.supabase.co")
CURRENT_KEY=$(get_value "SUPABASE_KEY" "your-supabase-anon-key")
CURRENT_SECRET=$(get_value "SECRET_KEY" "dev-secret-key-please-change-in-production")

# Prompt for new values
echo -e "${BLUE}Enter Supabase URL${NC} [default: $CURRENT_URL]:"
read -p "> " SUPABASE_URL
SUPABASE_URL=${SUPABASE_URL:-$CURRENT_URL}

echo -e "${BLUE}Enter Supabase Anonymous Key${NC} [default: $CURRENT_KEY]:"
read -p "> " SUPABASE_KEY
SUPABASE_KEY=${SUPABASE_KEY:-$CURRENT_KEY}

echo -e "${BLUE}Enter Secret Key for JWT encryption${NC} [default: $CURRENT_SECRET]:"
read -p "> " SECRET_KEY
SECRET_KEY=${SECRET_KEY:-$CURRENT_SECRET}

# Update .env file
cat > "$ENV_FILE" << EOF
# Supabase credentials (required for application to run)
SUPABASE_URL=$SUPABASE_URL
SUPABASE_KEY=$SUPABASE_KEY
SECRET_KEY=$SECRET_KEY

# Development settings
FLASK_APP=app.py
FLASK_ENV=development
LOG_LEVEL=DEBUG
LOG_TO_FILE=1
LOG_TO_CONSOLE=1
STRICT_SECRET_MODE=false
STRICT_ENV_MODE=false
EOF

echo -e "${GREEN}Environment configuration updated successfully!${NC}"
echo "Your .env file has been updated with the following values:"
echo -e "${BLUE}SUPABASE_URL:${NC} $SUPABASE_URL"
echo -e "${BLUE}SUPABASE_KEY:${NC} ${SUPABASE_KEY:0:8}..." # Show just first 8 chars for security
echo -e "${BLUE}SECRET_KEY:${NC} ${SECRET_KEY:0:8}..." # Show just first 8 chars for security
echo
echo "You're now ready to run the application with Podman:"
echo "cd .. && ./quickstart_podman.sh"