#!/bin/bash
# Status check script for CryoProtect containers

# Color definitions
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Container Status${NC}"
echo "================================="

# Check App container
echo -e "${BLUE}App Container (cryoprotect-app):${NC}"
if podman container exists "cryoprotect-app"; then
    if podman container inspect "cryoprotect-app" --format '{{.State.Running}}' | grep -q "true"; then
        echo -e "  Status: ${GREEN}Running${NC}"
        echo "  API URL: http://localhost:5001"
        echo "  Health Endpoint: http://localhost:5001/health"
    else
        echo -e "  Status: ${YELLOW}Stopped${NC}"
    fi
else
    echo -e "  Status: ${RED}Not Created${NC}"
fi

# Check RDKit container
echo -e "${BLUE}RDKit Container (cryoprotect-rdkit):${NC}"
if podman container exists "cryoprotect-rdkit"; then
    if podman container inspect "cryoprotect-rdkit" --format '{{.State.Running}}' | grep -q "true"; then
        echo -e "  Status: ${GREEN}Running${NC}"
        echo "  API URL: http://localhost:5002"
        echo "  Health Endpoint: http://localhost:5002/health"
        
        # Check RDKit status
        echo "  RDKit Status:"
        curl -s http://localhost:5002/health | jq || echo "  Could not get RDKit status"
    else
        echo -e "  Status: ${YELLOW}Stopped${NC}"
    fi
else
    echo -e "  Status: ${RED}Not Created${NC}"
fi

# Check container communication
echo -e "${BLUE}Container Communication:${NC}"
# Check if app can reach RDKit
APP_CAN_REACH_RDKIT=false
if podman container exists "cryoprotect-app" && podman container inspect "cryoprotect-app" --format '{{.State.Running}}' | grep -q "true"; then
    if curl -s http://localhost:5001/rdkit/check | jq -r '.rdkit_service.available' | grep -q "true"; then
        APP_CAN_REACH_RDKIT=true
        echo -e "  App → RDKit: ${GREEN}OK${NC}"
    else
        echo -e "  App → RDKit: ${RED}Failed${NC}"
    fi
else
    echo -e "  App → RDKit: ${YELLOW}App not running${NC}"
fi

# Print overall status
echo
echo -e "${BLUE}Overall Status:${NC}"
if podman container exists "cryoprotect-app" &&    podman container exists "cryoprotect-rdkit" &&    podman container inspect "cryoprotect-app" --format '{{.State.Running}}' | grep -q "true" &&    podman container inspect "cryoprotect-rdkit" --format '{{.State.Running}}' | grep -q "true" &&    [ "$APP_CAN_REACH_RDKIT" = "true" ]; then
    echo -e "  ${GREEN}All containers running and communicating properly${NC}"
else
    echo -e "  ${RED}Issues detected with container setup${NC}"
fi
