#!/bin/bash
# Main script to handle the complete repository transfer process

# Set base variables
BASE_DIR="/home/mushu/Projects/CryoProtect"
OLD_ORG="blueprint-house"
NEW_ORG="mushu-dev"
REPO_NAME="CryoProtect"

# Function to display help information
show_help() {
    echo "Repository Transfer Tool"
    echo "========================"
    echo "This script automates the process of transferring a GitHub repository"
    echo "from an organization to a personal account, including all necessary updates."
    echo ""
    echo "Usage: $0 [options] [steps]"
    echo ""
    echo "Options:"
    echo "  -b, --base-dir DIR     Set the base directory (default: $BASE_DIR)"
    echo "  -o, --old-org NAME     Set the old organization name (default: $OLD_ORG)"
    echo "  -n, --new-org NAME     Set the new organization name (default: $NEW_ORG)"
    echo "  -r, --repo-name NAME   Set the repository name (default: $REPO_NAME)"
    echo "  -d, --dry-run          Don't make changes, just show what would be done"
    echo "  -h, --help             Display this help message"
    echo ""
    echo "Steps:"
    echo "  all                    Run all steps (default)"
    echo "  urls                   Update GitHub URLs in project files"
    echo "  workflows              Update GitHub Actions workflow files"
    echo "  vercel                 Update Vercel deployment configuration"
    echo "  heroku                 Update Heroku deployment configuration"
    echo ""
    echo "Example:"
    echo "  $0 --dry-run              # Dry run all steps"
    echo "  $0 urls workflows         # Only update URLs and workflows"
    echo "  $0                        # Run all steps with default settings"
}

# Process command line options
DRY_RUN=0
STEPS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--base-dir)
            BASE_DIR="$2"
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
        -r|--repo-name)
            REPO_NAME="$2"
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
        all|urls|workflows|vercel|heroku)
            STEPS+=("$1")
            shift
            ;;
        *)
            echo "Unknown option or step: $1"
            show_help
            exit 1
            ;;
    esac
done

# If no steps specified, run all
if [ ${#STEPS[@]} -eq 0 ]; then
    STEPS=("all")
fi

# Validate base directory
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Base directory '$BASE_DIR' does not exist or is not a directory"
    exit 1
fi

# Change to base directory
cd "$BASE_DIR" || exit 1

# Display the current configuration
echo "Repository Transfer Configuration:"
echo "--------------------------------"
echo "Base directory:     $BASE_DIR"
echo "Old organization:   $OLD_ORG"
echo "New organization:   $NEW_ORG"
echo "Repository name:    $REPO_NAME"
echo "Dry run:            $([ $DRY_RUN -eq 1 ] && echo "Yes" || echo "No")"
echo "Steps to run:       ${STEPS[*]}"
echo "--------------------------------"

# Ask for confirmation
read -p "Do you want to proceed with these settings? (y/N) " CONFIRM
if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
    echo "Operation cancelled"
    exit 0
fi

# Function to run a step if it's included
run_step() {
    local step_name="$1"
    local step_command="$2"
    
    if [[ " ${STEPS[*]} " =~ " all " ]] || [[ " ${STEPS[*]} " =~ " $step_name " ]]; then
        echo ""
        echo "Running step: $step_name"
        echo "--------------------------------"
        
        if [ $DRY_RUN -eq 1 ]; then
            # Add --dry-run flag if the script supports it
            if [[ "$step_command" == *".sh"* ]]; then
                step_command="$step_command --dry-run"
            fi
        fi
        
        # Execute the command
        eval "$step_command"
        
        # Check the result
        if [ $? -eq 0 ]; then
            echo "Step completed: $step_name"
        else
            echo "Step failed: $step_name"
            echo "Do you want to continue with the next steps? (y/N)"
            read -p " " CONTINUE
            if [[ ! "$CONTINUE" =~ ^[Yy]$ ]]; then
                echo "Process aborted"
                exit 1
            fi
        fi
    fi
}

# Run the appropriate steps
run_step "urls" "./scripts/update-github-urls.sh --base-path \"$BASE_DIR\" --old-org \"$OLD_ORG\" --new-org \"$NEW_ORG\""
run_step "workflows" "./scripts/update-github-workflows.sh --base-path \"$BASE_DIR\" --old-org \"$OLD_ORG\" --new-org \"$NEW_ORG\""
run_step "vercel" "./scripts/update-vercel-repository.sh"
run_step "heroku" "./scripts/update-heroku-repository.sh"

echo ""
echo "Repository transfer process completed"
echo "--------------------------------"

if [ $DRY_RUN -eq 1 ]; then
    echo "This was a dry run, no actual changes were made."
    echo "Run without --dry-run to apply changes."
fi

echo ""
echo "IMPORTANT: Manual Steps Required"
echo "--------------------------------"
echo "1. Transfer the repository on GitHub:"
echo "   - Go to: https://github.com/$OLD_ORG/$REPO_NAME/settings"
echo "   - Scroll to the 'Danger Zone' section and click 'Transfer'"
echo "   - Enter '$OLD_ORG/$REPO_NAME' and '$NEW_ORG' when prompted"
echo ""
echo "2. After the transfer is complete, verify the results:"
echo "   - Check that the repository is available at: https://github.com/$NEW_ORG/$REPO_NAME"
echo "   - Verify the GitHub Actions workflows are set up correctly"
echo "   - Test the Vercel, Heroku, and Fly.io deployments"
echo ""
echo "3. Update GitHub Secrets:"
echo "   - Transfer all secrets from the old repository to the new one"
echo ""
echo "For more details, see the REPOSITORY_TRANSFER_GUIDE.md file"