#!/bin/bash
# Test IPv4 connectivity to Supabase from Podman container

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Supabase Connectivity Test${NC}"
echo "This script will test IPv4 connectivity to your Supabase instance from a Podman container"
echo "==========================================================================="

# Check if the Supabase URL is set in .env
if [ -f ../.env ]; then
  source ../.env
  echo -e "${BLUE}Using Supabase URL from .env:${NC} $SUPABASE_URL"
else
  echo -e "${YELLOW}No .env file found. Using default Supabase URL.${NC}"
  SUPABASE_URL="https://your-project.supabase.co"
fi

# Extract the hostname from the URL
SUPABASE_HOST=$(echo "$SUPABASE_URL" | sed -e 's|^[^/]*//||' -e 's|/.*$||')
echo -e "${BLUE}Extracted Supabase hostname:${NC} $SUPABASE_HOST"

# Test IPv4 connectivity from Podman container
echo -e "\n${BLUE}Testing IPv4 DNS Resolution from Podman Container${NC}"
echo "Attempting to resolve $SUPABASE_HOST using IPv4 from container..."
podman run --rm --network=cryoprotect-net alpine sh -c "apk add --no-cache bind-tools && host -4 $SUPABASE_HOST" || \
  (echo -e "${RED}Failed to resolve IPv4 from container${NC}" && DNS_CONTAINER_IPV4_FAILED=1)

# Test HTTP connectivity from Podman container
echo -e "\n${BLUE}Testing HTTP Connectivity from Podman Container${NC}"
echo "Attempting to connect to $SUPABASE_URL from container..."
podman run --rm --network=cryoprotect-net alpine sh -c "apk add --no-cache curl && curl -s -I \"$SUPABASE_URL\" | head -n 1" || \
  (echo -e "${RED}Failed to connect from container${NC}" && HTTP_CONTAINER_FAILED=1)

# Summary
echo -e "\n${BLUE}Connectivity Test Summary${NC}"
echo "==========================================================================="

if [ -n "$DNS_CONTAINER_IPV4_FAILED" ]; then
  echo -e "${RED}✗ Container IPv4 DNS resolution:${NC} Failed"
else
  echo -e "${GREEN}✓ Container IPv4 DNS resolution:${NC} Passed"
fi

if [ -n "$HTTP_CONTAINER_FAILED" ]; then
  echo -e "${RED}✗ Container HTTP connectivity:${NC} Failed"
else
  echo -e "${GREEN}✓ Container HTTP connectivity:${NC} Passed"
fi

# Recommendations
echo -e "\n${BLUE}Recommendations${NC}"
echo "==========================================================================="

if [ -n "$DNS_CONTAINER_IPV4_FAILED" ]; then
  echo -e "${YELLOW}Container DNS Resolution Issues:${NC}"
  echo "  - Try using a custom DNS server in your Podman container:"
  echo "    podman run --dns 1.1.1.1 --dns 8.8.8.8 ..."
  echo "  - Or edit podman-compose.yml to add DNS servers:"
  echo "      dns:"
  echo "        - 1.1.1.1"
  echo "        - 8.8.8.8"
fi

if [ -n "$HTTP_CONTAINER_FAILED" ]; then
  echo -e "${YELLOW}Container HTTP Connectivity Issues:${NC}"
  echo "  - Try running containers with host networking:"
  echo "    podman run --network=host ..."
  echo "  - Or update podman-compose.yml to use host networking mode:"
  echo "      network_mode: host"
fi

if [ -z "$DNS_CONTAINER_IPV4_FAILED" ] && [ -z "$HTTP_CONTAINER_FAILED" ]; then
  echo -e "${GREEN}Your Podman container can successfully connect to Supabase using IPv4.${NC}"
  echo "You should be able to proceed with running the application."
fi