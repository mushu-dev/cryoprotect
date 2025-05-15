#!/bin/bash
# Script to update GitHub Actions workflow files to use the new repository location

# Set variables
BASE_PATH="/home/mushu/Projects/CryoProtect"
WORKFLOWS_DIR="$BASE_PATH/.github/workflows"
OLD_ORG="blueprint-house"
NEW_ORG="mushu-dev"
LOG_FILE="$BASE_PATH/workflow_updates.log"

# Function to display help information
show_help() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -b, --base-path PATH   Set the base path (default: $BASE_PATH)"
    echo "  -o, --old-org NAME     Set the old organization name (default: $OLD_ORG)"
    echo "  -n, --new-org NAME     Set the new organization name (default: $NEW_ORG)"
    echo "  -d, --dry-run          Don't make changes, just show what would be done"
    echo "  -h, --help             Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 --dry-run              # Show changes without applying them"
    echo "  $0                        # Apply all changes"
}

# Process command line options
DRY_RUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--base-path)
            BASE_PATH="$2"
            WORKFLOWS_DIR="$BASE_PATH/.github/workflows"
            shift 2
            ;;
        -o|--old-org)
            OLD_ORG="$2"
            shift 2
            ;;
        -n|--new-org)
            NEW_ORG="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate paths
if [ ! -d "$WORKFLOWS_DIR" ]; then
    echo "Error: Workflows directory '$WORKFLOWS_DIR' does not exist or is not a directory"
    exit 1
fi

# Create log file
echo "Starting GitHub Actions workflow update $(date)" > "$LOG_FILE"
echo "From: $OLD_ORG to $NEW_ORG" >> "$LOG_FILE"
echo "Base path: $BASE_PATH" >> "$LOG_FILE"
echo "Workflows directory: $WORKFLOWS_DIR" >> "$LOG_FILE"
echo "Dry run: $([ $DRY_RUN -eq 1 ] && echo "Yes" || echo "No")" >> "$LOG_FILE"
echo "------------------------------------------------" >> "$LOG_FILE"

# Get list of workflow files
WORKFLOW_FILES=$(find "$WORKFLOWS_DIR" -type f -name "*.yml" -o -name "*.yaml")

if [ -z "$WORKFLOW_FILES" ]; then
    echo "No workflow files found in $WORKFLOWS_DIR"
    echo "No workflow files found in $WORKFLOWS_DIR" >> "$LOG_FILE"
    exit 0
fi

echo "Found $(echo "$WORKFLOW_FILES" | wc -l) workflow files"
echo "Found $(echo "$WORKFLOW_FILES" | wc -l) workflow files" >> "$LOG_FILE"

# Process each workflow file
for file in $WORKFLOW_FILES; do
    file_basename=$(basename "$file")
    echo ""
    echo "Processing: $file_basename"
    echo "Processing: $file" >> "$LOG_FILE"
    
    # Check if the file contains the old organization name
    if grep -q "$OLD_ORG" "$file"; then
        echo "Found references to $OLD_ORG in $file_basename"
        echo "Found references to $OLD_ORG in $file" >> "$LOG_FILE"
        
        # Show context of the matches
        grep -n "$OLD_ORG" "$file" >> "$LOG_FILE"
        
        # Apply the change if not in dry run mode
        if [ $DRY_RUN -eq 0 ]; then
            # Use sed to replace the pattern
            sed -i "s|$OLD_ORG|$NEW_ORG|g" "$file"
            echo "Updated: $file" >> "$LOG_FILE"
            echo "  ✓ Updated $file_basename"
        else
            echo "  (Dry run) Would update $file_basename"
        fi
    else
        echo "No references to $OLD_ORG found in $file_basename"
        echo "No references to $OLD_ORG found in $file" >> "$LOG_FILE"
    fi
    
    # Check for other GitHub-related settings
    if grep -q "github.com" "$file"; then
        echo "Found GitHub URL references in $file_basename:"
        grep -n "github.com" "$file" | grep -v "$NEW_ORG" >> "$LOG_FILE"
        grep -n "github.com" "$file" | grep -v "$NEW_ORG" 
        echo "Check these manually to ensure they're updated correctly"
    fi
done

echo ""
echo "------------------------------------------------"
echo "GitHub Actions workflow update completed"
echo "See $LOG_FILE for details"

if [ $DRY_RUN -eq 1 ]; then
    echo ""
    echo "This was a dry run, no actual changes were made."
    echo "Run without --dry-run to apply changes."
fi

# Instructions for GitHub Secrets
echo ""
echo "Don't forget to update GitHub Secrets for the new repository:"
echo "1. Go to https://github.com/$OLD_ORG/CryoProtect/settings/secrets/actions"
echo "2. Copy all secrets"
echo "3. Add them to https://github.com/$NEW_ORG/CryoProtect/settings/secrets/actions"
echo ""
echo "Required secrets include:"
echo "- HEROKU_API_KEY"
echo "- HEROKU_APP_NAME"
echo "- HEROKU_EMAIL"
echo "- VERCEL_ORG_ID"
echo "- VERCEL_PROJECT_ID"
echo "- VERCEL_TOKEN"
echo "- FLY_API_TOKEN"
echo "- Any other secrets used in your workflows"