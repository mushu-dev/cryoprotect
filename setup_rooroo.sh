#!/bin/bash
# Setup script for RooRoo GitHub Issue Manager in CryoProtect project

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}Setting up RooRoo GitHub Issue Manager for CryoProtect...${NC}"

# Check if we're in the right directory
if [ ! -d ".github/rooroo" ]; then
    echo -e "${RED}Error: .github/rooroo directory not found. Are you in the right directory?${NC}"
    echo "Please run this script from the root of the CryoProtect project."
    exit 1
fi

# Check if GitHub CLI is installed
if ! command -v gh &> /dev/null; then
    echo -e "${YELLOW}GitHub CLI not found. Installing...${NC}"
    
    # Check the OS
    if command -v dnf &> /dev/null; then
        # Fedora
        sudo dnf install -y gh
    elif command -v apt &> /dev/null; then
        # Ubuntu/Debian
        sudo apt update
        sudo apt install -y gh
    elif command -v pacman &> /dev/null; then
        # Arch Linux
        sudo pacman -S --noconfirm github-cli
    else
        echo -e "${RED}Unable to install GitHub CLI automatically. Please install it manually:${NC}"
        echo "https://github.com/cli/cli#installation"
        exit 1
    fi
    
    echo -e "${GREEN}GitHub CLI installed successfully.${NC}"
else
    echo -e "${GREEN}GitHub CLI already installed.${NC}"
fi

# Make RooRoo scripts executable
echo "Making RooRoo scripts executable..."
chmod +x .github/rooroo/issue-manager.sh
chmod +x .github/rooroo/scripts/*.sh
echo -e "${GREEN}Scripts are now executable.${NC}"

# Check GitHub authentication status
if ! gh auth status &> /dev/null; then
    echo -e "${YELLOW}GitHub CLI is not authenticated. Please authenticate:${NC}"
    echo "You will need to authenticate with GitHub to use RooRoo."
    echo "Choose 'GitHub.com' and 'HTTPS' when prompted."
    echo "Press Enter to continue with authentication..."
    read
    
    gh auth login
    
    # Check if authentication was successful
    if ! gh auth status &> /dev/null; then
        echo -e "${RED}GitHub authentication failed.${NC}"
        exit 1
    fi
else
    echo -e "${GREEN}GitHub CLI is already authenticated.${NC}"
fi

# Create RooRoo configuration directory
echo "Creating RooRoo configuration directory..."
mkdir -p .github/rooroo/config

# Create a basic configuration file
cat > .github/rooroo/config/config.yml << 'EOF'
repository: $(basename $(pwd))
base_branch: master
task_manager_id: cryoprotect-fedora-migration

# Labels to use
labels:
  - fedora-migration
  - selinux
  - database
  - docker
  - environment
  - testing
  - documentation
  - p0-blocker
  - p1-critical
  - p2-important
  - p3-enhancement

# Issue templates
templates:
  default:
    title: "[Fedora] {{title}}"
    body: |
      ## Objective
      {{description}}
      
      ## Steps
      1. {{steps}}
      
      ## Verification
      - [ ] Works on Fedora 38+
      - [ ] SELinux compatibility verified
      - [ ] Documentation updated
      
      ## Dependencies
      {{dependencies}}
EOF

echo -e "${GREEN}Configuration file created at .github/rooroo/config/config.yml${NC}"

# Create a test issue
echo "Do you want to create a test issue to verify RooRoo is working? (y/n)"
read -r CREATE_TEST_ISSUE

if [ "$CREATE_TEST_ISSUE" = "y" ] || [ "$CREATE_TEST_ISSUE" = "Y" ]; then
    echo "Creating a test issue..."
    .github/rooroo/issue-manager.sh create-issue \
      --title "Test Issue from Fedora Linux" \
      --body "This is a test issue created from Fedora Linux to verify RooRoo GitHub Issue Manager is working properly." \
      --labels "test,fedora-verification"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Test issue created successfully!${NC}"
    else
        echo -e "${RED}Failed to create test issue. Please check the error message above.${NC}"
    fi
fi

echo -e "${GREEN}RooRoo GitHub Issue Manager setup complete!${NC}"
echo 
echo "To use RooRoo, run commands like:"
echo ".github/rooroo/issue-manager.sh create-issue --title \"Issue Title\" --body \"Issue Body\" --labels \"label1,label2\""
echo
echo "For more information, see the documentation at:"
echo "https://github.com/rswaminathan/rooroo-github/blob/main/README.md"