#!/bin/bash
# Script to update GitHub URLs from blueprint-house to mushu-dev

# Set the base path
BASE_PATH="/home/mushu/Projects/CryoProtect"
OLD_ORG="blueprint-house"
NEW_ORG="mushu-dev"
LOG_FILE="$BASE_PATH/url_updates.log"

# Function to display help information
show_help() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -b, --base-path PATH   Set the base path to search (default: $BASE_PATH)"
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

# Validate base path
if [ ! -d "$BASE_PATH" ]; then
    echo "Error: Base path '$BASE_PATH' does not exist or is not a directory"
    exit 1
fi

# Create log file
echo "Starting URL update $(date)" > "$LOG_FILE"
echo "From: $OLD_ORG to $NEW_ORG" >> "$LOG_FILE"
echo "Base path: $BASE_PATH" >> "$LOG_FILE"
echo "Dry run: $([ $DRY_RUN -eq 1 ] && echo "Yes" || echo "No")" >> "$LOG_FILE"
echo "------------------------------------------------" >> "$LOG_FILE"

# Function to update files
update_files() {
    local search_pattern="$1"
    local replace_pattern="$2"
    local file_patterns="$3"
    
    echo "Searching for '$search_pattern' in files matching: $file_patterns"
    echo "Searching for '$search_pattern' in files matching: $file_patterns" >> "$LOG_FILE"
    
    # Use find to get all files matching the patterns
    found_files=$(find "$BASE_PATH" -type f -not -path "*/node_modules/*" -not -path "*/\.git/*" -not -path "*/\.vercel/*" -not -path "*/venv/*" -not -path "*/\.cache/*" -not -path "*/\__pycache__/*" $file_patterns)
    
    # Process each found file
    for file in $found_files; do
        # Check if the file contains the search pattern
        if grep -q "$search_pattern" "$file"; then
            echo "Found in: $file"
            echo "Found in: $file" >> "$LOG_FILE"
            
            # Show context of the matches
            grep -n "$search_pattern" "$file" >> "$LOG_FILE"
            
            # Apply the change if not in dry run mode
            if [ $DRY_RUN -eq 0 ]; then
                # Use sed to replace the pattern
                sed -i "s|$search_pattern|$replace_pattern|g" "$file"
                echo "Updated: $file" >> "$LOG_FILE"
                echo "  ✓ Updated"
            else
                echo "  (Dry run) Would update"
            fi
        fi
    done
}

# Update GitHub URLs in various file types
echo "Updating GitHub URLs from $OLD_ORG to $NEW_ORG..."

# Markdown files
update_files "$OLD_ORG" "$NEW_ORG" "-name \"*.md\""

# YAML files
update_files "$OLD_ORG" "$NEW_ORG" "-name \"*.yml\" -o -name \"*.yaml\""

# JavaScript/TypeScript files
update_files "$OLD_ORG" "$NEW_ORG" "-name \"*.js\" -o -name \"*.ts\" -o -name \"*.jsx\" -o -name \"*.tsx\""

# JSON files
update_files "$OLD_ORG" "$NEW_ORG" "-name \"*.json\""

# Python files
update_files "$OLD_ORG" "$NEW_ORG" "-name \"*.py\""

# Shell scripts
update_files "$OLD_ORG" "$NEW_ORG" "-name \"*.sh\" -o -name \"*.bash\""

# HTML files
update_files "$OLD_ORG" "$NEW_ORG" "-name \"*.html\" -o -name \"*.htm\""

# Update git remote URL
if [ $DRY_RUN -eq 0 ]; then
    echo "Updating Git remote URL..."
    cd "$BASE_PATH" || exit
    git remote set-url origin "https://github.com/$NEW_ORG/CryoProtect.git"
    echo "Git remote URL updated"
else
    echo "(Dry run) Would update Git remote URL to: https://github.com/$NEW_ORG/CryoProtect.git"
fi

echo "------------------------------------------------"
echo "URL update completed"
echo "See $LOG_FILE for details"

if [ $DRY_RUN -eq 1 ]; then
    echo ""
    echo "This was a dry run, no actual changes were made."
    echo "Run without --dry-run to apply changes."
fi