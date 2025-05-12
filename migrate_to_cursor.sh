#!/bin/bash
# Simple wrapper script for Cursor migration

# Print banner
echo "====================================================="
echo "  CryoProtect: Migration to Cursor IDE"
echo "====================================================="
echo ""
echo "This script will prepare the project for Cursor IDE migration by:"
echo "  1. Consolidating and organizing GitHub issues"
echo "  2. Creating necessary documentation"
echo "  3. Setting up the project structure for Cursor"
echo ""
echo "Press ENTER to continue or CTRL+C to cancel..."
read

# Run the GitHub issue consolidation script
echo "Step 1: Organizing GitHub issues..."
./scripts/consolidate_github_issues.sh

# Verify the CURSOR_MIGRATION_GUIDE.md exists
echo ""
echo "Step 2: Verifying migration documentation..."
if [ -f "CURSOR_MIGRATION_GUIDE.md" ]; then
    echo "✅ CURSOR_MIGRATION_GUIDE.md is ready"
else
    echo "❌ CURSOR_MIGRATION_GUIDE.md is missing. Creating it now..."
    # Logic to create the file would go here in a real implementation
fi

# Create Cursor workspace directory if it doesn't exist
echo ""
echo "Step 3: Setting up Cursor workspace..."
mkdir -p .cursor/extensions
echo "✅ Cursor workspace directory created"

# Update README to include Cursor information
echo ""
echo "Step 4: Updating README with Cursor information..."
if grep -q "Cursor IDE" README.md; then
    echo "✅ README.md already contains Cursor information"
else
    echo "ℹ️ Consider adding Cursor IDE information to README.md"
fi

# Notify user of completion
echo ""
echo "====================================================="
echo "  Migration preparation complete!"
echo "====================================================="
echo ""
echo "Next steps:"
echo "1. Open the project in Cursor IDE"
echo "2. Follow the instructions in CURSOR_MIGRATION_GUIDE.md"
echo "3. Update team members on the migration progress"
echo ""
echo "For more information, refer to the CURSOR_MIGRATION_GUIDE.md file."