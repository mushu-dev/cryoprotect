#!/bin/bash
# Script to list all CryoProtect containers and their purposes

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== CryoProtect Containers ====${NC}"
echo

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo -e "${RED}Error: Podman is not installed.${NC}"
    echo "Please install Podman first using:"
    echo "  sudo dnf install podman"
    exit 1
fi

# Get a list of all containers
containers=$(podman ps -a --format "{{.Names}}")

if [ -z "$containers" ]; then
    echo -e "${YELLOW}No containers found.${NC}"
    echo "You can create the RDKit container with:"
    echo "  ./quick_conda_container.sh"
    exit 0
fi

# Print a table header
printf "${CYAN}%-25s %-15s %-40s${NC}\n" "CONTAINER NAME" "STATUS" "PURPOSE"
printf "${CYAN}%-25s %-15s %-40s${NC}\n" "-------------" "------" "-------"

# Check each container
for container in $containers; do
    # Get container status
    status=$(podman container inspect -f '{{.State.Status}}' "$container")
    
    # Highlight our special RDKit container
    if [[ "$container" == "CryoProtect-RDKit-Conda" ]]; then
        name_color=$GREEN
        purpose="RDKit conda environment for molecular properties"
    elif [[ "$container" == *"RDKit"* ]]; then
        name_color=$GREEN
        purpose="Container with RDKit support"
    else
        name_color=$MAGENTA
        # Get container purpose from labels if available
        purpose=$(podman container inspect -f '{{.Config.Labels.purpose}}' "$container" 2>/dev/null || echo "Unknown")
        
        # If purpose is empty or null, try to guess based on image
        if [[ -z "$purpose" || "$purpose" == "<no value>" || "$purpose" == "null" ]]; then
            image=$(podman container inspect -f '{{.Config.Image}}' "$container")
            if [[ "$image" == *"miniconda"* || "$image" == *"conda"* ]]; then
                purpose="Conda environment container"
            elif [[ "$image" == *"slim"* ]]; then
                purpose="Slim/minimal container"
            else
                purpose="General purpose container"
            fi
        fi
    fi
    
    # Color status based on running/exited
    if [[ "$status" == "running" ]]; then
        status_color=$GREEN
    else
        status_color=$YELLOW
    fi
    
    # Print container info
    printf "${name_color}%-25s ${status_color}%-15s ${NC}%-40s\n" "$container" "$status" "$purpose"
done

echo
echo -e "${BLUE}=== Usage Instructions ===${NC}"
echo
echo -e "${GREEN}For RDKit and molecular property calculation:${NC}"
echo -e "  Use container: ${GREEN}CryoProtect-RDKit-Conda${NC}"
echo
echo "Run the unified ChEMBL import with RDKit:"
echo "  ./run_unified_import_in_container.sh"
echo
echo "Run any script with RDKit support:"
echo "  ./run_with_rdkit.sh your_script.py"
echo 
echo "Access the RDKit container directly:"
echo "  podman exec -it CryoProtect-RDKit-Conda bash -c 'source ~/.bashrc && conda activate cryoprotect && bash'"