#!/bin/bash
# make-scripts-executable.sh - Make all shell scripts executable
#
# This script makes all shell scripts in the scripts directory executable.
# It should be run after checking out the repository.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Making shell scripts executable..."

# Make all .sh files executable
find "$SCRIPT_DIR" -name "*.sh" -type f -exec chmod +x {} \;

echo "Done!"