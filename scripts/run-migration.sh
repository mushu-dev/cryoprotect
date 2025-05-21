#!/bin/bash

# Supabase to Convex Migration Script
# This script runs the migration with environment variables

set -e # Exit on error

# Colors for better output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Default values
CHECKPOINT_DIR="./migration_checkpoints"
LOG_FILE="./migration.log"

# Print header
echo -e "${BLUE}=======================================${NC}"
echo -e "${BLUE}  Supabase to Convex Migration Tool    ${NC}"
echo -e "${BLUE}=======================================${NC}"
echo ""

# Check if required environment variables are set
if [ -z "$SUPABASE_URL" ]; then
  echo -e "${YELLOW}SUPABASE_URL is not set. Please enter your Supabase URL:${NC}"
  read -r SUPABASE_URL
  export SUPABASE_URL
fi

if [ -z "$SUPABASE_KEY" ]; then
  echo -e "${YELLOW}SUPABASE_KEY is not set. Please enter your Supabase service_role key:${NC}"
  read -r SUPABASE_KEY
  export SUPABASE_KEY
fi

if [ -z "$CONVEX_URL" ]; then
  echo -e "${YELLOW}CONVEX_URL is not set. Please enter your Convex deployment URL:${NC}"
  read -r CONVEX_URL
  export CONVEX_URL
fi

# Create directories if they don't exist
mkdir -p "$CHECKPOINT_DIR"
mkdir -p "$(dirname "$LOG_FILE")"

# Run migration
echo -e "${BLUE}Starting migration from Supabase to Convex...${NC}"
echo -e "${YELLOW}This may take a while depending on the amount of data.${NC}"
echo -e "${YELLOW}You can check progress in the log file: ${LOG_FILE}${NC}"
echo ""

# Run the migration script
node ./scripts/migrate-to-convex.js

# Check if migration was successful
if [ $? -eq 0 ]; then
  echo -e "${GREEN}Migration completed successfully!${NC}"
  echo -e "${GREEN}The log file is available at: ${LOG_FILE}${NC}"
else
  echo -e "${RED}Migration failed. Please check the log file for details: ${LOG_FILE}${NC}"
  exit 1
fi