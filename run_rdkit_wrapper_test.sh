#!/bin/bash
# Test RDKit wrapper in both host and container environment

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Testing RDKit Wrapper Implementation${NC}"
echo "========================================="

# First test in host environment
echo -e "${BLUE}Testing in host environment:${NC}"
python3 rdkit_wrapper.py

# Make the wrapper executable
chmod +x rdkit_wrapper.py

# Now test in container environment if available
CONTAINER_NAME="CryoProtect-RDKit-Conda"

if podman container exists $CONTAINER_NAME && podman container inspect $CONTAINER_NAME --format '{{.State.Running}}' | grep -q "true"; then
    echo -e "\n${BLUE}Testing in container environment:${NC}"
    
    # Copy the wrapper and mock implementation to the container
    podman cp rdkit_wrapper.py $CONTAINER_NAME:/app/
    podman cp mock_rdkit_formula.py $CONTAINER_NAME:/app/
    
    # Run the test in the container
    podman exec -it $CONTAINER_NAME bash -c "cd /app && python3 rdkit_wrapper.py"
    
    echo -e "\n${GREEN}Tests completed in both environments${NC}"
else
    echo -e "\n${YELLOW}Container $CONTAINER_NAME not running, skipping container test${NC}"
    echo "Run './quick_conda_container.sh' first to set up the container"
    echo -e "${GREEN}Host environment test completed${NC}"
fi

# Print summary
echo -e "\n${BLUE}Integration Strategy:${NC}"
echo "1. Use the rdkit_wrapper.py module for all RDKit functionality"
echo "2. This provides a unified interface that works with both real RDKit and mock implementation"
echo "3. For production, deploy using the container with RDKit properly installed"
echo "4. See RDKIT_INTEGRATION_GUIDE.md for full implementation details"

# Make the script executable
chmod +x $(basename $0)