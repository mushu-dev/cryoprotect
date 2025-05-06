#!/bin/bash
# Test script for RooRoo GitHub Issue Manager on Fedora Linux

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

ROOROO_PATH=".github/rooroo/issue-manager.sh"

# Function to print section headers
print_section() {
    echo -e "\n${BLUE}======== $1 ========${NC}\n"
}

# Function to check if a test passes
check_result() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✅ PASS: $1${NC}"
        return 0
    else
        echo -e "${RED}❌ FAIL: $1${NC}"
        return 1
    fi
}

print_section "RooRoo GitHub Issue Manager Fedora Test Suite"

echo "This test suite will verify that RooRoo GitHub Issue Manager"
echo "is working correctly on your Fedora Linux system."
echo
echo -e "${YELLOW}Note: This will create actual GitHub issues in your repository.${NC}"
echo "These can be deleted after testing."
echo

read -p "Enter the repository owner/name (e.g., yourusername/CryoProtect): " REPO
if [ -z "$REPO" ]; then
    echo -e "${RED}Repository name is required. Exiting.${NC}"
    exit 1
fi

# Test 1: Basic Issue Creation
print_section "Test 1: Basic Issue Creation"
echo "Creating a simple issue..."

$ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Test Issue from Fedora" \
  --body "This is a test issue created from Fedora Linux to verify the RooRoo GitHub Issue Manager." \
  --labels "test,fedora-verification"

check_result "Basic issue creation"

# Test 2: Multi-line Content Handling
print_section "Test 2: Multi-line Content Handling"
echo "Creating an issue with multi-line content..."

$ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Multi-line Test from Fedora" \
  --body "This issue tests multi-line content handling on Fedora Linux.

Line 1: Testing line breaks
Line 2: Testing paragraph formatting
* Testing bullet points
* Testing nested content

## Testing markdown headers
\`\`\`bash
# Testing code blocks
echo 'Hello from Fedora'
\`\`\`

Verify formatting is preserved." \
  --labels "test,formatting-test"

check_result "Multi-line content handling"

# Test 3: Label Management
print_section "Test 3: Label Management"
echo "Creating fedora-specific labels..."

$ROOROO_PATH create-labels \
  --repo "$REPO" \
  --labels "fedora,selinux,database,docker,environment,testing"

check_result "Label creation"

echo "Creating an issue with multiple labels..."

$ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Label Test from Fedora" \
  --body "Testing label application from Fedora Linux" \
  --labels "fedora,selinux,test"

check_result "Multiple label application"

# Test 4: Issue Updates
print_section "Test 4: Issue Updates"
echo "Creating an issue and then updating it..."

ISSUE_NUMBER=$($ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Update Test from Fedora" \
  --body "Initial content" \
  --labels "test" \
  --output-number)

check_result "Issue creation for update test"

echo "Issue number: $ISSUE_NUMBER"
echo "Updating the issue..."

$ROOROO_PATH update-issue \
  --repo "$REPO" \
  --issue-number "$ISSUE_NUMBER" \
  --body "Updated content from Fedora Linux

Testing if updates preserve formatting and structure.
- Item 1
- Item 2

This update was made through RooRoo on Fedora."

check_result "Issue update"

# Test 5: SELinux Compatibility
print_section "Test 5: SELinux Compatibility"
echo "Checking SELinux status..."

if command -v getenforce &> /dev/null; then
    SELINUX_STATUS=$(getenforce)
    echo "SELinux status: $SELINUX_STATUS"
    
    if [ "$SELINUX_STATUS" = "Enforcing" ]; then
        echo "Creating issue with SELinux enforcing..."
        
        $ROOROO_PATH create-issue \
          --repo "$REPO" \
          --title "SELinux Enforcing Test" \
          --body "This issue was created while SELinux was in enforcing mode on Fedora." \
          --labels "test,selinux"
        
        check_result "Issue creation with SELinux enforcing"
    else
        echo -e "${YELLOW}SELinux is not in enforcing mode. Skipping this test.${NC}"
        echo "To fully test SELinux compatibility, run with: sudo setenforce 1"
    fi
else
    echo -e "${YELLOW}SELinux tools not found. Skipping this test.${NC}"
fi

# Test 6: Error Recovery
print_section "Test 6: Error Recovery"
echo "Testing error recovery by attempting to update a non-existent issue..."

$ROOROO_PATH update-issue \
  --repo "$REPO" \
  --issue-number "999999999" \
  --body "This update should fail gracefully" \
  || true

echo -e "${YELLOW}The error above is expected. Checking if RooRoo handles errors gracefully...${NC}"

# Test if RooRoo is still functioning after an error
$ROOROO_PATH create-issue \
  --repo "$REPO" \
  --title "Post-Error Test" \
  --body "This issue tests if RooRoo recovers after an error." \
  --labels "test"

check_result "Recovery after error"

# Final summary
print_section "Test Results Summary"

echo -e "${GREEN}All tests completed.${NC}"
echo 
echo "If all tests passed, RooRoo is working correctly on your Fedora system."
echo "You can now use RooRoo for managing GitHub issues for the CryoProtect project."
echo
echo "To use RooRoo, run commands like:"
echo "$ROOROO_PATH create-issue --repo \"$REPO\" --title \"Issue Title\" --body \"Issue Body\" --labels \"label1,label2\""
echo
echo -e "${YELLOW}Don't forget to delete the test issues if they're no longer needed.${NC}"
echo "You can delete them from the GitHub web interface."